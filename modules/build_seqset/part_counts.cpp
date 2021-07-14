#include "modules/build_seqset/part_counts.h"

namespace build_seqset {

std::pair<size_t, size_t> part_counts::seq_to_index_range(dna_slice seq) {
  size_t start_idx = 0;
  size_t last_idx = 0;

  for (size_t i = 0; i < m_bases; ++i) {
    start_idx <<= 2;
    last_idx <<= 2;
    if (i >= seq.size()) {
      last_idx += 3;
    } else {
      int b = int(seq[i]);
      start_idx += b;
      last_idx += b;
    }
  }

  size_t limit_idx = last_idx + 1;

  return std::make_pair(start_idx, limit_idx);
}

part_counts::part_counts(unsigned n_bases)
    : m_bases(n_bases), m_counts(track_alloc("part_counts")) {
  CHECK_LT(n_bases, sizeof(size_t) * 4);
  CHECK_LE(m_bases, seq_repository::k_inline_bases);
  m_counts.resize(1ULL << (2 * n_bases), 0);
}

std::string part_counts::display_histo() const {
  if (m_counts.empty()) {
    return "(no counts)";
  }
  std::stringstream out;

  tracked_vector<size_t> dist(m_counts);
  std::sort(dist.begin(), dist.end());

  size_t tot = 0;
  for (size_t val : dist) {
    tot += val;
  }

  out << "Total: " << tot << " Count: " << dist.size()
      << " Avg: " << (tot + dist.size() / 2) / dist.size();

  for (size_t pct : {0, 1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99, 100}) {
    out << " " << pct << "%: ";

    size_t pos = ((dist.size() - 1) * pct + 50) / 100;
    CHECK_LT(pos, dist.size());
    out << dist[pos];
  }

  return out.str();
}

}  // namespace build_seqset
