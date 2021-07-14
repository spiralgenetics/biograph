#include "modules/bio_base/seqset_flat.h"
#include "modules/io/parallel.h"
#include "modules/io/spiral_file_mem.h"

#include <deque>

const product_version k_seqset_flat_version{"1.0.0"};

seqset_flat_builder::seqset_flat_builder(const seqset* the_seqset)
    : m_seqset(the_seqset),
      m_whole_seqset_range(the_seqset->ctx_begin()),
      m_seqs("", track_alloc("seqset_flat:seqs")) {}

uint64_t seqset_flat_builder::add_sequence(const dna_slice& seq_part1,
                                           const dna_slice& seq_part2) {
  std::lock_guard<std::mutex> l(m_seq_mu);

  size_t orig_offset = m_seqs_offset;
  for (const dna_slice& seq : {seq_part1, seq_part2}) {
    for (dna_base base : seq) {
      m_queued_bases <<= 2;
      m_queued_bases |= int(base);
      m_seqs_offset++;
      if ((m_seqs_offset & 3) == 0) {
        m_seqs.write((const char*)&m_queued_bases, 1);
      }
    }
  }

  m_flat_seqs++;

  return orig_offset;
}

namespace {

// Enable extra debugging checks when producing a seqset_flat.  Causes a
// significant slowdown.
constexpr bool k_debug_flat = false;

constexpr uint64_t k_no_entry = std::numeric_limits<uint64_t>::max();

}  // namespace

void seqset_flat_builder::add_entry(uint64_t seqset_entry_id,
                                    uint64_t base_offset, uint64_t rel_offset,
                                    bool rc, const dna_sequence& seq_buffer) {
  uint64_t abs_offset = base_offset + rel_offset;
  CHECK_LE(abs_offset, k_max_flat_offset);

  if (rc) {
    DCHECK_GE(rel_offset + 1, m_seqset->entry_size(seqset_entry_id));
  }

  if (k_debug_flat) {
    dna_sequence sub;
    if (rc) {
      CHECK_LT(rel_offset, seq_buffer.size());
      sub = seq_buffer
                .subseq(rel_offset + 1 - m_seqset->entry_size(seqset_entry_id),
                        m_seqset->entry_size(seqset_entry_id))
                .rev_comp();
    } else {
      CHECK_LE(rel_offset + m_seqset->entry_size(seqset_entry_id),
               seq_buffer.size());
      sub =
          seq_buffer.subseq(rel_offset, m_seqset->entry_size(seqset_entry_id));
    }
    dna_sequence entry_seq = m_seqset->ctx_entry(seqset_entry_id).sequence();
    CHECK_EQ(sub.as_string(), entry_seq.as_string());
  }

  uint64_t offset_with_flags = (abs_offset << k_flag_bits);
  if (rc) {
    offset_with_flags |= k_rc_flag;
  }

  m_entries->set(seqset_entry_id, offset_with_flags);
}

void seqset_flat_builder::finalize() {
  if (m_seqs_offset & 3) {
    while (m_seqs_offset & 3) {
      m_queued_bases <<= 2;
      m_seqs_offset++;
    }
    m_seqs.write((const char*)&m_queued_bases, 1);
  }
  SPLOG(
      "%ld flattened sequence buffers written.  Total bases: %ld. "
      "Total size: %ld MB",
      m_flat_seqs, m_seqs_offset, m_seqs_offset / (4 * 1024 * 1024));
}

bool seqset_flat_builder::claim_entry(uint64_t seqset_entry_id) {
  if (m_claimed_entries->at(seqset_entry_id).safe_increment()) {
    // Already claimed.
    return false;
  }
  return true;
}

// State for tracing through a seqset, accumulating a flat sequence.
struct seqset_flat_builder::trace_state {
  // Sequence gathered so far.
  dna_sequence seq;
  // Number of bases still needed to be traced in order to fully
  // flatten sequences that we've claimed.
  unsigned fwd_needed = 0;

  // Current seqset entry id being traced with pop_front
  uint64_t pop_entry_id = k_no_entry;
  // Current seqset entry range being traced with push_front.
  seqset_range rc_range;

  // Entry ids we've claimed in the forward direction (from pop_front)
  std::vector<uint64_t> fwd_entries;
  // Entry ids we've claimed in the reverse complement direction (from
  // push_front)
  std::vector<uint64_t> rc_entries;

  // The entry we've seen when tracing in the reverse complement direction that
  // we can use to switch directions.
  uint64_t first_rc_entry = k_no_entry;
  unsigned first_rc_offset = 0;
};

// Traces one base in the forward direction.
void seqset_flat_builder::trace(trace_state& state) {
  if (claim_entry(state.pop_entry_id)) {
    state.fwd_entries.push_back(state.pop_entry_id);
    state.fwd_needed = std::max(
        state.fwd_needed, unsigned(m_seqset->entry_size(state.pop_entry_id)));
  } else {
    if (state.fwd_needed == 0) {
      // Exit early if we weren't able to claim anything at all.
      return;
    }
    state.fwd_entries.push_back(k_no_entry);
  }

  dna_base base = m_seqset->entry_get_base(state.pop_entry_id);
  uint64_t new_pop_entry_id = m_seqset->entry_pop_front(state.pop_entry_id);
  state.seq.push_back(base);
  state.fwd_needed--;

  state.rc_range = state.rc_range.push_front_drop(base.complement());
  CHECK(state.rc_range.valid());
  uint64_t rc_entry_id = state.rc_range.begin();
  if (state.rc_range.size() == state.seq.size()) {
    state.first_rc_entry = rc_entry_id;
    state.first_rc_offset = state.seq.size();
  }
  if (state.rc_range.size() == m_seqset->entry_size(rc_entry_id) &&
      claim_entry(rc_entry_id)) {
    state.rc_entries.push_back(rc_entry_id);
  } else {
    state.rc_entries.push_back(k_no_entry);
  }

  state.pop_entry_id = new_pop_entry_id;
}

uint64_t seqset_flat_builder::process_entry(uint64_t entry_id) {
  // State for tracing to the right, i.e. with
  // entry_pop_front(entry_id).
  trace_state right_trace;

  // State for tracing to the left, i.e. with
  // entry_pop_front(entry_id's reverse complement).
  trace_state left_trace;

  right_trace.rc_range = m_whole_seqset_range;
  left_trace.rc_range = m_whole_seqset_range;

  // Trace right until we've flattened all the sequences needed for
  // the entries we've claimed.
  right_trace.pop_entry_id = entry_id;
  right_trace.fwd_needed = 0;

  do {
    trace(right_trace);
  } while (right_trace.fwd_needed);
  if (right_trace.fwd_entries.empty()) {
    // Unable to claim anything.
    return 0;
  }

  // Trace left from the the beginning of our trace to the right.
  CHECK_NE(right_trace.first_rc_entry, k_no_entry);
  unsigned left_right_shared = right_trace.first_rc_offset;
  left_trace.pop_entry_id = right_trace.first_rc_entry;
  left_trace.fwd_needed = left_right_shared;

  CHECK_GT(left_trace.fwd_needed, 0);
  do {
    trace(left_trace);
  } while (left_trace.fwd_needed);

  // Construct a sequence buffer with the bases we've found from tracing left
  // and right.
  size_t right_offset = left_trace.seq.size() - left_right_shared;
  size_t left_offset = left_trace.seq.size() - 1;

  dna_slice left_seq_buffer = dna_slice(left_trace.seq).rev_comp();
  dna_slice right_seq_buffer =
      dna_slice(right_trace.seq)
          .subseq(left_right_shared,
                  right_trace.seq.size() - left_right_shared);
  size_t base_offset = add_sequence(left_seq_buffer, right_seq_buffer);

  dna_sequence debug_seq_buffer;
  if (k_debug_flat) {
    debug_seq_buffer =
        dna_sequence(left_seq_buffer.begin(), left_seq_buffer.end());
    debug_seq_buffer +=
        dna_sequence(right_seq_buffer.begin(), right_seq_buffer.end());
  }

  // Record pointers to our sequence buffer in all the entries we've
  // claimed.
  size_t num_processed = 0;
  for (bool trace_dir_left : {false, true}) {
    // trace_dir_left is true if we're storing sequences generated from tracing
    // left (through left_trace), and false if we're storing sequences
    const trace_state& state = trace_dir_left ? left_trace : right_trace;
    size_t offset_for_dir = trace_dir_left ? left_offset : right_offset;

    CHECK_EQ(state.fwd_entries.size(), state.rc_entries.size());
    CHECK_EQ(state.seq.size(), state.fwd_entries.size());
    for (bool trace_rc : {false, true}) {
      // trace_rc is true if we're storing entries found by calling push_front,
      // or false if we're storing entries found by calling pop_front.
      const std::vector<uint64_t>& entries =
          trace_rc ? state.rc_entries : state.fwd_entries;
      // If we're tracing left and we get a rc sequence by calling push_front,
      // the direction of the flattened sequence in our buffer will be to the
      // right.
      bool is_rc = trace_dir_left ? !trace_rc : trace_rc;

      for (unsigned i = 0; i < entries.size(); ++i) {
        uint64_t entry_id = entries[i];
        if (entry_id != k_no_entry) {
          num_processed++;
          add_entry(
              entry_id, base_offset,
              trace_dir_left ? (offset_for_dir - i) : (offset_for_dir + i),
              is_rc, debug_seq_buffer);
        }
      }
    }
  }
  return num_processed;
}

void seqset_flat_builder::build(const spiral_file_create_state& state,
                                progress_handler_t progress) {
  state.set_version("seqset_flat", k_seqset_flat_version);

  seqset_flat_metadata metadata;
  metadata.seqset_uuid = m_seqset->uuid();
  state.create_json<seqset_flat_metadata>("seqset_flat.json", metadata);

  m_entries.reset(new mutable_packed_varint_vector(
      state.create_subpart("entry_index",
                           state.options().with_delayed_write(true)),
      m_seqset->size(),
      (k_max_flat_offset << k_flag_bits) | ((1ULL << k_flag_bits) - 1)));
  m_claimed_entries.reset(new mutable_packed_vector<unsigned, 1>(
      m_seqset->size(), "seqset_flat_builder_claimed"));

  uint64_t tot_processed(0);
  std::mutex stats_mutex;
  parallel_for(0, m_seqset->size(), [this, &progress, &tot_processed,
                                     &stats_mutex](size_t start, size_t limit) {
    uint64_t chunk_entries_processed = 0;
    for (size_t idx = start; idx != limit; idx++) {
      chunk_entries_processed += process_entry(idx);

      if (chunk_entries_processed > (limit - start)) {
        std::lock_guard<std::mutex> lock(stats_mutex);
        // Checkpoint progress if we've traced through a lot.
        tot_processed += chunk_entries_processed;
        progress(tot_processed * 1.0 / m_seqset->size());
        chunk_entries_processed = 0;
      }
    }
  });
  progress(1.0);

  finalize();

  membuf seqs(new borrowed_membuf(m_seqs.buffer(), m_seqs.size()));
  state.create_membuf("sequence_data", seqs);
}

seqset_flat::seqset_flat(const spiral_file_open_state& state,
                         const seqset* the_seqset)
    : m_seqset(the_seqset) {
  state.enforce_max_version("seqset_flat", k_seqset_flat_version);

  m_metadata = state.open_json<seqset_flat_metadata>("seqset_flat.json");
  CHECK_EQ(m_metadata.seqset_uuid, the_seqset->uuid());
  m_seqs = state.open_membuf("sequence_data");
  m_entries.reset(new packed_varint_vector(state.open_subpart("entry_index")));
}
