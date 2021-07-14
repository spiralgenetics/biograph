#include <parallel/algorithm>
#include <endian.h>

#include "base/base.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_mapred/make_readmap.h"
#include "modules/io/binary_find.h"
#include "modules/io/file_io.h"
#include "modules/io/parallel.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/mapred/manifest_parallel.h"

void make_readmap::do_make(const std::string& readmap_file_path, const seqset_file& the_seqset_file,
                           manifest corrected_reads, bool is_paired, unsigned max_read_len,
                           progress_handler_t progress) {
  make_readmap maker(&the_seqset_file);
  spiral_file_create_mmap c(readmap_file_path);
  maker.create_from_reads(corrected_reads, is_paired, max_read_len, c.create(),
                          progress);
}

void make_readmap::migrate(const seqset_file& old_seqset,
                           const readmap& old_readmap,
                           const seqset_file& new_seqset,
                           const std::string& new_readmap_path,
                           bool throw_on_read_not_in_new_seqset) {
  make_readmap maker(&new_seqset);
  spiral_file_create_mmap c(new_readmap_path);
  maker.create_from_migrate(old_seqset, old_readmap, c.create(),
                            throw_on_read_not_in_new_seqset);
}

void make_readmap::upgrade(const readmap& old_readmap,
                                     const seqset_file& the_seqset_file,
                                     const std::string& new_readmap_path,
          const std::function<dna_sequence(uint64_t /* seqset_id */, unsigned)>&
          lookup_seq,
                                     progress_handler_t progress) {
  make_readmap maker(&the_seqset_file);
  spiral_file_create_mmap c(new_readmap_path);
  maker.create_from_upgrade(old_readmap, the_seqset_file, c.create(), lookup_seq, progress);
}

void make_readmap::fast_migrate(const readmap& old_readmap,
                                const seqset_mergemap& mergemap,
                                const spiral_file_create_state& new_readmap,
                                progress_handler_t progress) {
  make_readmap maker(nullptr);
  maker.create_from_fast_migrate(old_readmap, mergemap, new_readmap, progress);
}

void make_readmap::import_reads_from(
	 manifest corrected_reads_manifest
     , bool is_paired
     , progress_handler_t progress)
{
	SPLOG("make_readmap::do_make> Creating readmap, is_paired = %u", unsigned(is_paired));

    progress(0);

	manifest_reader corrected_reads_reader(corrected_reads_manifest);
	kv_reader corrected_reads_kvs(corrected_reads_reader);
	std::string read_id;
	corrected_reads read_pair;

	SPLOG("Constructing mate table builder.");
	size_t mate_loop_table_size = (is_paired ? 4 : 2) * corrected_reads_manifest.get_num_records();
    CHECK(m_mate_loop_table.empty());
    m_mate_loop_table.resize(mate_loop_table_size);
	SPLOG("Allocated %lu entries at %lu bytes each, total memory = %lu",
          m_mate_loop_table.size(), sizeof(mate_loop_table_entry)
		, m_mate_loop_table.size() * sizeof(mate_loop_table_entry)
	);

	parallel_mate_loop_table_builder the_table_builder(
		m_mate_loop_table
		, m_seqset
		, corrected_reads_manifest.count_file_infos()
		, corrected_reads_manifest.get_num_records()
		, is_paired
	);
	SPLOG("Starting mate loop table build.");
    progress(0.2);
	auto functor
		 = manifest_parallelize<parallel_mate_loop_table_builder, std::string, corrected_reads>
        (corrected_reads_manifest, the_table_builder, subprogress(progress, 0.2, 0.4));
	auto counts = functor.get_counts();
	SPLOG("Mate loop table has %lu entries, %lu paired and %lu unpaired", m_mate_loop_table.size(), counts.first, counts.second);
    progress(0.4);

	parallel_sort_in_place(m_mate_loop_table.begin(), m_mate_loop_table.end());
	SPLOG("Mate pair table is sorted.");
    progress(0.6);

	auto last_non_empty_ri
		= std::find_if_not(
			m_mate_loop_table.rbegin()
			, m_mate_loop_table.rend()
			, [](const mate_loop_table_entry& entry) { return entry.entry_id == k_no_loop_entry; }
		);
    progress(0.8);
	if (last_non_empty_ri == m_mate_loop_table.rend()) {
        m_mate_loop_table.clear();
        return;
    }
	m_mate_loop_table.erase(last_non_empty_ri.base(), m_mate_loop_table.end());
	SPLOG("Mate loop table empty entries dropped. Length = %lu.", m_mate_loop_table.size());
    progress(1);
}

std::pair<uint64_t, uint64_t> make_readmap::parallel_mate_loop_table_builder::get_counts() const
{
	return std::make_pair(
		std::accumulate(m_paired_counts.begin(), m_paired_counts.end(), 0ULL)
		, std::accumulate(m_unpaired_counts.begin(), m_unpaired_counts.end(), 0ULL)
	);
}

void make_readmap::parallel_mate_loop_table_builder::operator()(
    const std::string& read_id, const corrected_reads& read_pair,
    size_t file_info_id, size_t record_id) {
  dna_slice sequence = read_pair[0].corrected;
  dna_slice mate_sequence;
  if (read_pair.size() == 2 && m_is_paired) {
    mate_sequence = read_pair[1].corrected;
    m_paired_counts[file_info_id] += m_is_paired ? 4 : 2;
  } else if (read_pair.size() == 1) {
    m_unpaired_counts[file_info_id] += m_is_paired ? 4 : 2;
  } else {
    throw io_exception(
        boost::format("Unexpected read pairing found for read \"%1%\": %2% "
                      "reads were found in a \"pair.\" m_is_paired = %3%") %
        read_id % read_pair.size() % m_is_paired);
  }

  auto get_entry = [&](dna_slice read_sequence,
                       const char* desc) -> uint64_t {
    // Cautious readmap lookup is significantly slower, but will find
    // internal errors where the seqset and readmap don't match.
#if NDEBUG
    constexpr bool k_cautious_readmap_lookup = false;
#else
    constexpr bool k_cautious_readmap_lookup = true;
#endif
    if (k_cautious_readmap_lookup) {
      seqset_range r = m_seqset->find(read_sequence);
      if (!r.valid()) {
        throw io_exception(printstring(
            "Read record ID %lu, \"%s\" (%s) was not found in seqset.", record_id,
            desc, read_sequence.as_string().c_str()));
      }
      CHECK_EQ(r.begin(), m_seqset->find_existing_unique(read_sequence, 20));
      return r.begin();
    } else {
      return m_seqset->find_existing_unique(read_sequence, 20);
    }
  };

  unsigned read_len = sequence.size();
  uint64_t entry_id = get_entry(sequence, "forward");
  uint64_t rc_entry_id =
      get_entry(sequence.rev_comp(), "reverse");

  unsigned mate_read_len = mate_sequence.size();
  if (mate_read_len) {
    uint64_t mate_entry_id = get_entry(mate_sequence, "mate forward");
    uint64_t mate_rc_entry_id =
        get_entry(mate_sequence.rev_comp(), "mate reverse");

    if (sequence > mate_sequence) {
      // Canonicalize so our generated readmap is deterministic.
      std::swap(mate_entry_id, entry_id);
      std::swap(mate_read_len, read_len);
      std::swap(mate_rc_entry_id, rc_entry_id);
    }

    add_read_to_mate_loop_table(LOOP_START, record_id, entry_id, read_len, rc_entry_id, 0);
    add_read_to_mate_loop_table(RC, record_id, rc_entry_id, read_len,
                                mate_entry_id, mate_read_len);
    add_read_to_mate_loop_table(MATE, record_id, mate_entry_id, mate_read_len,
                                mate_rc_entry_id, 0);
    add_read_to_mate_loop_table(MATE_RC, record_id, mate_rc_entry_id,
                                mate_read_len, k_no_loop_entry, 0);
  } else {
    add_read_to_mate_loop_table(LOOP_START, record_id, entry_id, read_len, rc_entry_id, 0);
    add_read_to_mate_loop_table(RC, record_id, rc_entry_id, read_len, k_no_loop_entry, 0);
  }
}

void make_readmap::parallel_mate_loop_table_builder::
    add_read_to_mate_loop_table(loop_entry_type type, size_t record_id,
                                uint64_t entry_id, unsigned len,
                                uint64_t loop_entry_id, unsigned mate_len) {
  if (!m_is_paired) {
    CHECK_LT(int(type), 2);
  }
  if (type != RC) {
    CHECK_EQ(0, mate_len);
  }
  if (type == MATE_RC) {
    CHECK_EQ(k_no_loop_entry, loop_entry_id);
  } else if (type == LOOP_START || type == MATE) {
    CHECK_NE(k_no_loop_entry, loop_entry_id);
  }
  m_mate_loop_table[calc_mate_loop_table_id(record_id, m_is_paired) +
                    int(type)] =
      mate_loop_table_entry(type, entry_id, len, loop_entry_id, mate_len);
}

void make_readmap::create_common(const spiral_file_create_state& state,
                                 const std::string& seqset_uuid,
                                 size_t seqset_size,
                                 size_t num_reads, unsigned max_read_len) {
  state.set_version("readmap", k_readmap_version);

  readmap_metadata metadata;
  metadata.seqset_uuid = seqset_uuid;
  state.create_json<readmap_metadata>("readmap.json", metadata);

  m_sparse_multi.reset(new sparse_multi_builder(
      state.create_subpart("read_ids"), seqset_size, num_reads));

  m_read_lengths = make_unique<mutable_packed_varbit_vector>(
      state.create_subpart("read_lengths"), num_reads, max_read_len);
}

void make_readmap::create_from_reads(manifest corrected_reads_manifest, bool is_paired,
                                     unsigned max_read_len, const spiral_file_create_state& state,
                                     progress_handler_t progress) {
  import_reads_from(corrected_reads_manifest, is_paired, subprogress(progress, 0, 0.7));

  create_common(state, m_seqset->uuid(), m_seqset->size(), m_mate_loop_table.size(), max_read_len);

  for (uint64_t idx = 0; idx < m_mate_loop_table.size(); ++idx) {
    if ((idx & 0xFFFF) == 0) {
      subprogress(progress, 0.7, 0.75)(idx * 1. / m_mate_loop_table.size());
    }
    const auto& mate_loop_row = m_mate_loop_table[idx];
    m_sparse_multi->add(mate_loop_row.entry_id);
  }
  m_sparse_multi->finalize();

  progress(0.8);

  SPLOG("Filling read lengths");
  for (size_t read_id = 0; read_id < m_mate_loop_table.size(); ++read_id) {
    m_read_lengths->set(read_id, m_mate_loop_table[read_id].read_length);
  }
  m_pairing_data_present = true;
  if (m_pairing_data_present) {
    SPLOG("Processing pairing data");
    m_mate_loop_ptr.reset(new mutable_packed_varbit_vector(
        state.create_subpart("mate_loop_ptr"), m_mate_loop_table.size(), m_mate_loop_table.size()));
    m_is_forward.reset(new mutable_packed_vector<unsigned, 1>(state.create_subpart("is_forward"),
                                                              m_mate_loop_table.size()));

    // Finds the first matching read of the given specification.
    auto find_first_of = [&](loop_entry_type type, uint64_t entry_id,
                             unsigned read_length) -> uint64_t {
      auto try_it = std::lower_bound(m_mate_loop_table.begin(), m_mate_loop_table.end(),
                                     mate_loop_table_entry(type, entry_id, read_length, 0, 0));
      size_t try_idx = try_it - m_mate_loop_table.begin();
      return try_idx;
    };

    // First pass: search for RCs and mates and save the first mathcing entry.
    SPLOG("Filling mate loop entries in parallel");
    parallel_for(
        0, m_mate_loop_table.size(),
        [&](size_t start, size_t limit) {
          for (size_t idx = start; idx != limit; ++idx) {
            auto& loop_row = m_mate_loop_table[idx];
            switch (loop_row.type) {
              case LOOP_START:
                // save RC
                m_is_forward->at(idx) = 1;
                m_mate_loop_ptr->set(
                    idx, find_first_of(RC, loop_row.loop_entry_id, loop_row.read_length));
                break;
              case RC:
                if (loop_row.loop_entry_id != k_no_loop_entry) {
                  // Save mate.
                  m_mate_loop_ptr->set(
                      idx, find_first_of(MATE, loop_row.loop_entry_id, loop_row.mate_read_length));
                }
                break;
              case MATE:
                // save RC
                m_is_forward->at(idx) = 1;
                m_mate_loop_ptr->set(
                    idx, find_first_of(MATE_RC, loop_row.loop_entry_id, loop_row.read_length));
                break;
              case MATE_RC:
                // Fill in later when doing claims.
                break;
            }
          }
        },
        subprogress(progress, 0.8, 0.9));

    // Second pass, nonparallelizable so we can stay deterministic.  Claim entries and complete the
    // link from MATE_RC->LOOP_START

    SPLOG("Linking mate loops");
    subprogress fill_progress(progress, 0.9, 1);
    mutable_packed_vector<unsigned, 1> claimed(m_mate_loop_table.size(), "make_readmap:claimed");
    // Given the result from find_first_of, claims the next unclaimed read.  To ensure determinism,
    // this part should be run without parallelism.
    auto claim_next = [&](uint64_t try_idx, loop_entry_type type, uint64_t entry_id,
                          unsigned read_length) -> uint64_t {
      CHECK_LT(try_idx, m_mate_loop_table.size());
      {
        const auto& try_loop_row = m_mate_loop_table[try_idx];
        DCHECK_EQ(try_loop_row.entry_id, entry_id);
        DCHECK_EQ(try_loop_row.read_length, read_length);
        DCHECK_EQ(try_loop_row.type, type);
      }

      try_idx = claimed.claim_next_available(try_idx);
      CHECK_LT(try_idx, m_mate_loop_table.size());

      const auto& try_loop_row = m_mate_loop_table[try_idx];
      CHECK_EQ(try_loop_row.entry_id, entry_id);
      CHECK_EQ(try_loop_row.read_length, read_length);
      CHECK_EQ(try_loop_row.type, type);

      return try_idx;
    };

    for (size_t idx = 0; idx != m_mate_loop_table.size(); ++idx) {
      auto& loop_row = m_mate_loop_table[idx];
      if (loop_row.type != LOOP_START) {
        continue;
      }

      if ((idx & 0xFFFF) == 0) {
        fill_progress(idx * 1.0 / m_mate_loop_table.size());
      }

      uint64_t rc_idx =
          claim_next(m_mate_loop_ptr->get(idx), RC, loop_row.loop_entry_id, loop_row.read_length);
      m_mate_loop_ptr->set(idx, rc_idx);
      auto& rc_loop_row = m_mate_loop_table[rc_idx];

      if (rc_loop_row.loop_entry_id == k_no_loop_entry) {
        // No mate for this one; just point back to original.
        m_mate_loop_ptr->set(rc_idx, idx);
        continue;
      }

      uint64_t mate_idx = claim_next(m_mate_loop_ptr->get(rc_idx), MATE, rc_loop_row.loop_entry_id,
                                     rc_loop_row.mate_read_length);
      m_mate_loop_ptr->set(rc_idx, mate_idx);
      auto& mate_loop_row = m_mate_loop_table[mate_idx];

      uint64_t rc_mate_idx = claim_next(m_mate_loop_ptr->get(mate_idx), MATE_RC,
                                        mate_loop_row.loop_entry_id, rc_loop_row.mate_read_length);
      m_mate_loop_ptr->set(mate_idx, rc_mate_idx);
      // Save loop back to the beginning
      m_mate_loop_ptr->set(rc_mate_idx, idx);
    }

    SPLOG("Mate loop entries complete");
  }
  progress(1);
}

void make_readmap::create_from_migrate(const seqset_file& old_seqset_file,
                                       const readmap& old_readmap,
                                       const spiral_file_create_state& state,
                                       bool throw_on_read_not_in_new_seqset) {
  create_common(state, m_seqset->uuid(), m_seqset->size(),
                old_readmap.m_sparse_multi->dest_elem_count(), old_readmap.max_read_len());

  const auto& old_seqset = old_seqset_file.get_seqset();
  old_seqset.populate_pop_front_cache();

  size_t read_idx = 0;
  std::vector<uint8_t> read_lengths;
  read_lengths.resize(old_readmap.m_sparse_multi->dest_elem_count());
  for (const auto& entry : *old_readmap.m_sparse_multi) {
    uint64_t seqset_entry_id = entry.first;
    uint64_t read_id_range_start = entry.second.first;
    uint64_t read_id_range_end = entry.second.second;
    seqset_range old_seqset_entry_range = old_seqset.ctx_entry(seqset_entry_id);
    dna_sequence entry_seq = old_seqset_entry_range.sequence();
    seqset_range new_seqset_entry = m_seqset->find(entry_seq);
    if (throw_on_read_not_in_new_seqset && !new_seqset_entry.valid()) {
      throw io_exception(make_error_string(old_seqset_file, *m_seqset_file, entry_seq));
    }
    for (uint64_t read_id = read_id_range_start; read_id != read_id_range_end; read_id++) {
      m_sparse_multi->add(new_seqset_entry.begin());
      read_lengths[read_idx] = new_seqset_entry.size();
      ++read_idx;
    }
  }
  m_sparse_multi->finalize();

  for (size_t read_id = 0; read_id < read_lengths.size(); ++read_id) {
    m_read_lengths->set(read_id, read_lengths[read_id]);
  }
}

void make_readmap::create_from_upgrade(
    const readmap& old_readmap, const seqset_file& the_seqset_file,
    const spiral_file_create_state& new_readmap,
    const std::function<dna_sequence(uint64_t /* seqset_id */, unsigned)>& lookup_seq,
    progress_handler_t progress) {
  CHECK(old_readmap.has_pairing_data());
  const seqset* the_seqset = &the_seqset_file.get_seqset();
  create_common(new_readmap, old_readmap.metadata().seqset_uuid, the_seqset->size(),
                old_readmap.size(), old_readmap.max_read_len());

  m_pairing_data_present = true;

  for (size_t read_id = 0; read_id < old_readmap.size(); ++read_id) {
    m_read_lengths->set(read_id, old_readmap.get_readlength(read_id));
  }

  if (!old_readmap.has_mate_loop()) {
    const_cast<readmap&>(old_readmap)
        .enable_mate_loop(lookup_seq, subprogress(progress, 0.2, 0.7));
  }

  size_t total_added = 0;
  subprogress translate_progress(progress, 0.7, 0.9);
  for (const auto& entry : *old_readmap.m_sparse_multi) {
    for (uint64_t i = entry.second.first; i != entry.second.second; ++i) {
      CHECK_EQ(i, m_sparse_multi->add(entry.first));

      total_added++;
      if ((total_added & 65535) == 0) {
        translate_progress(total_added * 1.0 / old_readmap.size());
      }
    }
  }

  CHECK_EQ(total_added, old_readmap.size());

  m_sparse_multi->finalize();

  progress(0.95);
  m_mate_loop_ptr.reset(new mutable_packed_varbit_vector(
      new_readmap.create_subpart("mate_loop_ptr"), old_readmap.m_read_lengths->size(),
      old_readmap.m_read_lengths->size()));

  m_is_forward.reset(new mutable_packed_vector<unsigned, 1>(
      new_readmap.create_subpart("is_forward"), old_readmap.m_read_lengths->size()));

  progress(0.97);
  for (size_t i = 0; i < old_readmap.m_read_lengths->size(); ++i) {
    m_mate_loop_ptr->set(i, old_readmap.m_mate_loop_ptr->get(i));
    m_is_forward->at(i) = old_readmap.m_is_forward->at(i);
  }
}

void make_readmap::create_from_fast_migrate(
    const readmap& old_readmap, const seqset_mergemap& mergemap,
    const spiral_file_create_state& new_readmap,
    progress_handler_t progress) {
  CHECK_EQ(old_readmap.metadata().seqset_uuid, mergemap.metadata().orig_seqset_uuid);
  create_common(new_readmap, mergemap.metadata().merged_seqset_uuid,
                mergemap.get_bitcount().size(),
                old_readmap.m_sparse_multi->dest_elem_count(), old_readmap.max_read_len());

  for (size_t read_id = 0; read_id < old_readmap.size(); ++read_id) {
    m_read_lengths->set(read_id, old_readmap.get_readlength(read_id));
  }

  subprogress subp(progress, 0.2, 1);
  size_t total_added = 0;
  for (std::pair<uint64_t, std::pair<uint64_t, uint64_t>> entry :
       *old_readmap.m_sparse_multi) {
    uint64_t translated = mergemap.get_bitcount().find_count(entry.first);
    for (uint64_t i = entry.second.first; i != entry.second.second; ++i) {
      CHECK_EQ(i, m_sparse_multi->add(translated));

      total_added++;
      if ((total_added & 65535) == 0) {
        subp(total_added * 1.0 / old_readmap.size());
      }
    }
  }

  CHECK_EQ(total_added, old_readmap.size());

  m_sparse_multi->finalize();

  if (old_readmap.has_pairing_data()) {
    m_pairing_data_present = true;

    if (old_readmap.m_mate_loop_ptr) {
      CHECK_EQ(old_readmap.m_mate_loop_ptr->size(), old_readmap.m_read_lengths->size());
      m_mate_loop_ptr.reset(new mutable_packed_varbit_vector(
          new_readmap.create_subpart("mate_loop_ptr"),
          old_readmap.m_read_lengths->size(), old_readmap.m_read_lengths->size()));
    }
    if (old_readmap.m_mate_pair_ptr) {
      CHECK_EQ(old_readmap.m_mate_pair_ptr->size(), old_readmap.m_read_lengths->size());
      m_mate_pair_ptr.reset(new mutable_packed_vector<uint32_t, 32>(
          new_readmap.create_subpart("mate_pair_ptr"),
          old_readmap.m_read_lengths->size()));
    }
    CHECK_EQ(old_readmap.m_is_forward->size(), old_readmap.m_read_lengths->size());
    m_is_forward.reset(new mutable_packed_vector<unsigned, 1>(
        new_readmap.create_subpart("is_forward"),
        old_readmap.m_read_lengths->size()));

    for (size_t i = 0; i < old_readmap.m_read_lengths->size(); ++i) {
      if (old_readmap.m_mate_loop_ptr) {
        m_mate_loop_ptr->set(i, old_readmap.m_mate_loop_ptr->get(i));
      } else {
        m_mate_pair_ptr->at(i) = old_readmap.m_mate_pair_ptr->at(i);
      }
      m_is_forward->at(i) = old_readmap.m_is_forward->at(i);
    }
  }
}

std::string make_readmap::make_error_string(
	const seqset_file& old_seqset_file
	, const seqset_file& new_seqset_file
	, const dna_sequence&  entry_seq
)
{
	std::string error_string = "Readmap migration error while migrating "
		"readmap from source seqset \"";
	error_string += old_seqset_file.path();
	error_string += "\" to destination seqset \"";
	error_string += new_seqset_file.path();
	error_string += ". Sequence \"";
	error_string += entry_seq.as_string();
	error_string += "\" was present in the source seqset, but not in the destination."
		" Is the destination seqset a superset of the source seqset?";

	return error_string;
}
