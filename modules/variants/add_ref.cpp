#include "modules/variants/add_ref.h"
#include "modules/variants/scaffold.h"

namespace variants {

const char add_ref::k_add_ref_name[] = "ADD_REF";

add_ref::add_ref(const assemble_options& options, aoffset_t pad_size, bool whole_ref,
                 aoffset_t max_len, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)),
      m_cur_offset(-pad_size),
      m_pad_size(pad_size),
      m_max_len(max_len),
      m_options(options) {
  CHECK_GE(m_pad_size, 0);
  if (whole_ref) {
    m_padded_right_offset = std::numeric_limits<aoffset_t>::max() - m_pad_size - 1;
  }
}

void add_ref::output_ref(aoffset_t left_offset, aoffset_t right_offset) {
  CHECK_GT(right_offset, left_offset);

  if (left_offset < 0) {
    left_offset = 0;
  }

  if (m_max_len) {
    while (right_offset - left_offset > m_max_len) {
      output_ref_part(left_offset, left_offset + m_max_len);
      left_offset += m_max_len;
    }
  }

  if (right_offset <= left_offset) {
    return;
  }
  output_ref_part(left_offset, right_offset);
}

void add_ref::output_ref_part(aoffset_t left_offset, aoffset_t right_offset) {
  scaffold subs = m_options.scaffold->subscaffold(left_offset, right_offset - left_offset);
  // For most cases there will be only one extent in this scaffold
  for (const auto& ext : subs.extents()) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->matches_reference = true;
    a->tags.insert(k_add_ref_name);
    a->left_offset = left_offset + ext.offset;
    a->right_offset = left_offset + ext.offset + ext.sequence.size();
    a->seq = ext.sequence;
    sort_and_output(std::move(a));
  }
}

void add_ref::advance_to(aoffset_t left_offset) {
  aoffset_t padded_left_offset = left_offset - m_pad_size;
  while (!m_edge_offsets.empty() && *m_edge_offsets.begin() <= padded_left_offset) {
    aoffset_t next_edge = *m_edge_offsets.begin();
    if (next_edge > m_cur_offset) {
      output_ref(m_cur_offset, next_edge);
      m_cur_offset = next_edge;

      if (m_cur_offset > 0) {
        flush_sorted_to(m_cur_offset);
      }
    } else {
      CHECK_EQ(next_edge, m_cur_offset);
    }
    m_edge_offsets.erase(m_edge_offsets.begin());
  }

  if (!m_edge_offsets.empty()) {
    // No need to skip or emit padding right now.
    return;
  }

  if (padded_left_offset <= m_padded_right_offset) {
    // Haven't advanced far enough to need to emit padding.
    return;
  }

  if (m_cur_offset < m_padded_right_offset) {
    output_ref(m_cur_offset, m_padded_right_offset);
    m_cur_offset = m_padded_right_offset;
  }

  // Gap between m_padded_right_offset and padded_left_offset that we don't need to fill; skip it.
  CHECK_LT(m_cur_offset, padded_left_offset);
  m_cur_offset = padded_left_offset;

  if (m_cur_offset > 0) {
    flush_sorted_to(m_cur_offset);
  }
}

void add_ref::on_assembly(assembly_ptr a) {
  aoffset_t leftmost = min(a->left_offset, a->right_offset);
  advance_to(leftmost);

  if (m_cur_offset != leftmost) {
    m_edge_offsets.insert(leftmost);
  }
  if (a->right_offset) {
    m_edge_offsets.insert(a->right_offset);
    m_padded_right_offset = std::max<aoffset_t>(m_padded_right_offset, a->right_offset + m_pad_size);
  } else {
    m_padded_right_offset = std::max<aoffset_t>(m_padded_right_offset, a->left_offset + m_pad_size);
  }

  sort_and_output(std::move(a));
}

void add_ref::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_edge_offsets.empty());
}

}  // namespace variants
