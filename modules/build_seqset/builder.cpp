#include "modules/build_seqset/builder.h"
#include "modules/io/config.h"
#include "modules/io/parallel.h"
#include "modules/io/track_mem.h"

namespace build_seqset {

void builder::build_chunks(part_repo& entries, const std::string& pass_name,
                           bool keep_tmp, progress_handler_t progress) {
  SPLOG("Computing seqset parts");
  auto parts = entries.partitions(pass_name, true /* need pushed iterators */,
                                  not keep_tmp /* delete on close */);

  // Fixup for prefixes:
  // If we have a partition depth of 3, the front for "ACGT" might be "GA",
  // which would be in the "GAA" partition, not the "GAC" partition.
  for (unsigned i = 0; i + 1 < parts.size(); ++i) {
    auto& part = parts[i];
    auto& next = parts[i + 1];

    if (part.main && part.main->begin() != part.main->end()) {
      continue;
    }
    for (dna_base b : dna_bases()) {
      if (part.pushed[b].first != part.pushed[b].second) {
        std::swap(part.pushed[b], next.pushed[b]);
        std::swap(part.pushed_repositories[b], next.pushed_repositories[b]);
        CHECK(part.pushed[b].first == part.pushed[b].second);
      }
    }
  }

  std::mutex mu;

  dna_base_array<seq_repository::popped_iterator> whole_base_start,
      whole_base_end;

  std::map<unsigned, size_t> shared_histo;

  parallel_for(  //
      0, parts.size(),
      [&](size_t part_id) {
        const part_repo::partition_ref& part = parts[part_id];
        seq_repository::iterator cur_it = part.main->begin();
        seq_repository::iterator end_it = part.main->end();

        if (cur_it == end_it) {
          parts[part_id].reset();
          return;
        }

        dna_base_array<seq_repository::popped_iterator> base_cur_it;
        dna_base_array<seq_repository::popped_iterator> base_end_it;

        for (dna_base b : dna_bases()) {
          base_end_it[b] = part.pushed[b].second.pop_front();
          base_cur_it[b] = part.pushed[b].first.pop_front();
        }

        bool first_entry = true;
        auto prev_it = cur_it;

        std::unique_lock<std::mutex> l(mu);
        built_chunk& chunk = m_chunks[part.prefix];
        l.unlock();

        std::map<uint16_t, size_t> local_histo;
        unsigned local_max_read_len = 0;

        tracked_vector<uint16_t> sizes(track_alloc("builder:build_sizes_chunk"));
        tracked_vector<uint16_t> shared(track_alloc("builder:build_shared_chunk"));
        while (cur_it != end_it) {
          if (first_entry) {
            shared.push_back(0);
            first_entry = false;
          } else {
            unsigned shared_bases =
                cur_it->shared_prefix_length(*prev_it);
            shared.push_back(shared_bases);
            local_histo[shared_bases]++;
          }
          sizes.push_back(cur_it->size());
          local_max_read_len = std::max<unsigned>(local_max_read_len, cur_it->size());

          for (dna_base b : dna_bases()) {
            auto& bcur = base_cur_it[b];
            auto& bend = base_end_it[b];
            auto& prev = *chunk.has_prev[b];

            if (bcur == bend) {
              prev.push_back(false);
            } else {
              dna_compare_result cmp = bcur->compare_to(*cur_it);
              switch (cmp) {
                case dna_compare_result::FIRST_IS_LESS:
                  LOG(FATAL) << "Missing expansion?";
                case dna_compare_result::FIRST_IS_PREFIX:
                case dna_compare_result::EQUAL:
                  prev.push_back(true);
                  ++bcur;
                  break;
                case dna_compare_result::SECOND_IS_LESS:
                case dna_compare_result::SECOND_IS_PREFIX:
                  prev.push_back(false);
                  break;
              }
            }
          }

          prev_it = cur_it;
          ++cur_it;
        }

        for (dna_base b : dna_bases()) {
          CHECK(base_cur_it[b] == base_end_it[b]) << "base: " << b << " "
                                                  << base_cur_it[b]->sequence();
        }

        parts[part_id].reset();

        CHECK_EQ(sizes.size(), shared.size());
        if (sizes.empty()) {
          return;
        }
        chunk.sizes.emplace(sizes.size(), local_max_read_len, "build_seqset:sizes_chunk");
        chunk.shared.emplace(shared.size(), local_max_read_len, "build_seqset:shared_chunk");
        for (size_t i = 0; i != sizes.size(); ++i) {
          chunk.sizes->set(i, sizes[i]);
        }
        for (size_t i = 0; i != shared.size(); ++i) {
          chunk.shared->set(i, shared[i]);
        }

        l.lock();
        for (const auto& h : local_histo) {
          shared_histo[h.first] += h.second;
        }
        m_max_read_len = std::max<unsigned>(m_max_read_len, local_max_read_len);
      },
      progress);

  SPLOG("Maximum entry size: %d", m_max_read_len);
  SPLOG("Shared prefix histogram:");
  size_t tot_histo_entries = 0;
  for (const auto& h : shared_histo) {
    tot_histo_entries += h.second;
  }
  std::string line;
  size_t tot_histo_entries_so_far = 0;
  for (const auto& h : shared_histo) {
    tot_histo_entries_so_far += h.second;
    line += printstring(" %5d: %10ld (+%6.2f=%6.2f)", h.first, h.second,
                        h.second * 100. / tot_histo_entries,
                        tot_histo_entries_so_far * 100. / tot_histo_entries);
    if (line.size() > 100) {
      SPLOG("%s", line.c_str());
      line.clear();
    }
  }
  if (!line.empty()) {
    SPLOG("%s", line.c_str());
    line.clear();
  }
}

namespace {

void fill_prev(const tracked_vector<bool>& in, bitcount& out, size_t size, size_t out_offset) {
  auto in_it = in.begin();

  size_t out_end = out_offset + size;

  constexpr size_t k_pad_bits = 8 * sizeof(uint64_t);
  size_t in_offset = 0;
  while (out_offset < out_end && in_offset < k_pad_bits) {
    if (*in_it) {
      out.set(out_offset, true);
    }
    ++in_it;
    ++out_offset;
    ++in_offset;
  }

  // If we're not too close to either edge of the block, we can avoid
  // the atomic ops because we don't have to worry about interfering
  // with other blocks.
  while (out_offset + k_pad_bits < out_end) {
    if (*in_it) {
      out.set_unlocked(out_offset, true);
    }
    ++in_it;
    ++out_offset;
  }

  while (in_it != in.end()) {
    if (*in_it) {
      out.set(out_offset, true);
    }
    ++in_it;
    ++out_offset;
  }
  CHECK_EQ(out_offset, out_end);
}

}  // namespace

std::unique_ptr<seqset> builder::make_seqset(
    const spiral_file_create_state& state, progress_handler_t progress) {
  SPLOG("Calculating seqset size");

  std::vector<size_t> tot_offsets;
  std::vector<dna_sequence> chunk_parts;

  dna_sequence prev_chunk;
  size_t tot_size = 0;
  for (auto& chunk_out : m_chunks) {
    if (!chunk_out.second.sizes) {
      continue;
    }
    const dna_sequence& prefix = chunk_out.first;
    auto& chunk = chunk_out.second;
    tot_offsets.push_back(tot_size);
    chunk_parts.push_back(prefix);
    tot_size += chunk.sizes->size();
    chunk.shared->set(0, prefix.shared_prefix_length(prev_chunk));
    prev_chunk = prefix;
  }

  SPLOG("%ld total seqset entries; initializing seqset", tot_size);

  std::unique_ptr<seqset> result(new seqset(state, tot_size, m_max_read_len));
  result->init();

  SPLOG("Filling in seqset");

  parallel_for(  //
      0, chunk_parts.size(),
      [&](size_t idx) {
        const builder::built_chunk& chunk = m_chunks[chunk_parts[idx]];
        size_t chunk_offset = tot_offsets[idx];
        size_t tot_pos = chunk_offset;
        size_t chunk_size = chunk.sizes->size();

        for (size_t chunk_pos = 0; chunk_pos != chunk_size; chunk_pos++) {
          result->set_shared(tot_pos + chunk_pos, chunk.shared->get(chunk_pos));
        }
        for (size_t chunk_pos = 0; chunk_pos != chunk_size; chunk_pos++) {
          result->set_entry_size(tot_pos + chunk_pos, chunk.sizes->get(chunk_pos));
        }
        for (dna_base b : dna_bases()) {
          const auto& prev = *chunk.has_prev[b];
          auto& seqset_prev = result->mutable_prev(b);
          fill_prev(prev, seqset_prev, chunk_size, tot_pos);
        }
        tot_pos += chunk_size;
      },
      subprogress(progress, 0.2, 1));
  m_chunks.clear();

  SPLOG("Finalizing seqset");
  result->finalize();
  return result;
}

}  // namespace build_seqset
