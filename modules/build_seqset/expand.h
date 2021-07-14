#pragma once

#include "modules/build_seqset/part_repo.h"
#include "modules/build_seqset/part_counts.h"
#include "modules/io/progress.h"

namespace build_seqset {

class expander {
 public:
  expander(part_repo& entries, bool keep_tmp)
      : m_entries(entries), m_keep_tmp(keep_tmp) {}

  // Returns number of entries deduplicated
  size_t sort_and_dedup(const std::string& already_sorted_pass,
                        const std::string& new_entries_pass,
                        const std::string& result_sorted_pass,
                        const std::string& result_expanded_pass,
                        unsigned stride, unsigned count,
                        progress_handler_t progress = null_progress_handler);

  // Returns number of expansions done.
  //
  // Stride and count control what expansions are emitted once a
  // needed expansion is detected.
  //
  // If a popped-front entry is needed for "ABCDEFG" (e.g. "BCDEFG"),
  // here are what will be generated with various settings:
  //
  // stride=1 count=1:
  // BCDEFG
  // stride=2 count=1:
  // BCDEFG
  // stride=1 count=2:
  // BCDEFG
  // CDEFG
  // stride=2 count=2:
  // BCDEFG
  // DEFG
  // stride=2 count=255:
  // BCDEFG
  // DEFG
  // FG
  size_t expand(const std::string& input_pass, const std::string& expanded_pass,
                unsigned stride, unsigned count,
                progress_handler_t progress = null_progress_handler);

 private:
  part_repo& m_entries;
  bool m_keep_tmp;

  std::unique_ptr<part_counts> m_part_counts;
};

}  // namespace build_seqset
