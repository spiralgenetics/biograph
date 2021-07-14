#include "modules/variants/assemble_testutil.h"

namespace variants {

void PrintTo(const assembly& as, std::ostream* os) {
  *os << "\n";
  as.output_offsets(*os);
  as.output_other_info(*os);
  *os << ": ";

  aoffset_t left_anchor_len = as.left_anchor_len;
  aoffset_t right_anchor_len = as.right_anchor_len;
  if (left_anchor_len > aoffset_t(as.seq.size())) {
    left_anchor_len = as.seq.size();
  }
  if (right_anchor_len > aoffset_t(as.seq.size())) {
    right_anchor_len = as.seq.size();
  }
  aoffset_t right_anchor_start = aoffset_t(as.seq.size()) - right_anchor_len;

  dna_sequence left, main, right;
  if (int(as.seq.size()) >= (left_anchor_len + right_anchor_len)) {
    left = as.seq.subseq(0, left_anchor_len);
    right = as.seq.subseq(right_anchor_start, right_anchor_len);
    main = as.seq.subseq(left_anchor_len, right_anchor_start - left_anchor_len);
  } else {
    *os << "(ol) ";
    left = as.seq.subseq(0, right_anchor_start);
    right = as.seq.subseq(left_anchor_len, aoffset_t(as.seq.size()) - left_anchor_len);
    main = as.seq.subseq(right_anchor_start,
                         left_anchor_len + right_anchor_len - aoffset_t(as.seq.size()));
  }

  *os << left << " " << main << " " << right << " id=" << as.assembly_id;
  int num_ids = 0;
  for (const auto id : as.merged_assembly_ids) {
    *os << "," << id;
    if (++num_ids > 2) {
      *os << "...";
      break;
    }
  }
  *os << " score=" << as.score;
  if (as.strand_count) {
    *os << " strand_count=" << as.strand_count;
  }
  if (as.edge_coverage) {
    *os << " edge_coverage(" << *as.edge_coverage << ")";
  }

  if (!as.sub_assemblies.empty()) {
    *os << ", " << as.sub_assemblies.size() << " subassemblies for phase_ids=(";
    bool first = true;
    for (const auto& pid : as.phase_ids) {
      if (!first) {
        *os << ",";
      }
      first = false;
      *os << pid;
    }
    *os << "):\n";
    size_t n = 0;
    for (const auto& suba : as.sub_assemblies) {
      *os << "  #" << n << ": ";
      ++n;
      (*suba)->output_offsets(*os);
      *os << ": " << (*suba)->seq << "\n";
    }
  }
  *os << "\n";
}

std::vector<dna_sequence> reads_for_seq(dna_sequence seq, unsigned read_length,
                                        unsigned read_distance) {
  std::vector<dna_sequence> reads;
  CHECK_LE(read_length, seq.size());

  unsigned i = 0;
  while (i + read_length <= seq.size()) {
    reads.push_back(seq.subseq(i, read_length));
    i += read_distance;
  }

  // In case seq.size() - read_length is not an exact multiple of
  // read_distance, fill in another read at the end.
  reads.push_back(seq.subseq(seq.size() - read_length, read_length));
  return reads;
}

}  // namespace variants
