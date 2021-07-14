#include "modules/variants/scaffold.h"

namespace variants {
void scaffold::add(aoffset_t offset, dna_slice seq) {
  if (!m_extents.empty()) {
    CHECK_GE(offset, m_end_pos);
  }
  extent new_part;
  new_part.offset = offset;
  new_part.sequence = seq;
  m_extents.push_back(new_part);

  aoffset_t new_end = new_part.offset + aoffset_t(seq.size());
  if (new_end > m_end_pos) {
    m_end_pos = new_end;
  }
}

aoffset_t scaffold::calc_end_pos() const {
  if (m_extents.empty()) {
    return 0;
  }
  return m_extents.back().offset + aoffset_t(m_extents.back().sequence.size());
}

dna_slice scaffold::save_storage(const dna_sequence& seq) {
  if (!m_seq_storage) {
    m_seq_storage = std::make_shared<std::list<dna_sequence>>();
  }
  return *m_seq_storage->insert(m_seq_storage->begin(), seq);
}

bool scaffold::is_simple() const { return m_extents.size() == 1 && m_extents.front().offset == 0; }

dna_slice scaffold::get_simple() const {
  CHECK(is_simple());
  return m_extents.front().sequence;
}

scaffold scaffold::subscaffold(aoffset_t start, aoffset_t len) const {
  aoffset_t limit = start + len;

  CHECK_GE(start, 0);
  CHECK_GE(len, 0);
  CHECK_GE(limit, 0);

  scaffold result;
  result.m_seq_storage = m_seq_storage;
  for (const auto& part : m_extents) {
    if (part.offset >= limit) {
      continue;
    }
    if (part.offset + aoffset_t(part.sequence.size()) <= start) {
      continue;
    }

    aoffset_t relative_start = start - part.offset;
    if (relative_start < 0) {
      relative_start = 0;
    }
    aoffset_t relative_limit = limit - part.offset;
    if (relative_limit > aoffset_t(part.sequence.size())) {
      relative_limit = aoffset_t(part.sequence.size());
    }

    extent new_part;
    new_part.offset = (part.offset + relative_start) - start;
    new_part.sequence = part.sequence.subseq(relative_start, relative_limit - relative_start);
    result.m_extents.push_back(new_part);
  }
  result.m_end_pos = len;
  CHECK_GE(result.m_end_pos, result.calc_end_pos());
  return result;
}

std::string scaffold::subscaffold_str(aoffset_t start, aoffset_t len) const {
  aoffset_t limit = start + len;

  CHECK_GE(len, 0);

  std::string result;

  if (limit < 0) {
    return std::string(len, 'N');
  }
  aoffset_t cur_offset = 0;

  if (start < 0) {
    result += std::string(-start, 'N');
    len -= (-start);
    start = 0;
  }

  CHECK_GE(len, 0);
  CHECK_EQ(limit, start + len);

  for (const auto& part : m_extents) {
    if (part.offset >= limit) {
      continue;
    }
    if (part.offset + aoffset_t(part.sequence.size()) <= start) {
      continue;
    }

    aoffset_t relative_start = start - part.offset;
    if (relative_start < 0) {
      relative_start = 0;
    }
    aoffset_t relative_limit = limit - part.offset;
    if (relative_limit > aoffset_t(part.sequence.size())) {
      relative_limit = aoffset_t(part.sequence.size());
    }

    aoffset_t new_part_offset = (part.offset + relative_start) - start;
    CHECK_GE(new_part_offset, cur_offset);
    result += std::string(new_part_offset - cur_offset, 'N');
    cur_offset = new_part_offset;

    result += part.sequence.subseq(relative_start, relative_limit - relative_start).as_string();
    cur_offset += relative_limit - relative_start;
  }
  CHECK_LE(cur_offset, len);
  result += std::string(len - cur_offset, 'N');
  return result;
}

void scaffold::set_end_pos(aoffset_t new_end_pos) {
  CHECK_GE(new_end_pos, m_end_pos);
  m_end_pos = new_end_pos;
}

void scaffold::print_to(std::ostream& os) const {
  os << "Scaffold, [0-" << m_end_pos << "), " << m_extents.size() << " extent";
  if (m_extents.size() != 1) {
    os << "s";
  }
  os << ":\n";
  for (const auto& e : m_extents) {
    os << " [" << e.offset << "-" << e.offset + e.sequence.size() << "): " << e.sequence << "\n";
  }
  os << "Scaffold end at " << m_end_pos << "\n";
}

unsigned scaffold::shared_prefix_length(dna_slice seq) const {
  if (empty()) {
    return 0;
  }
  const auto& first = m_extents[0];
  if (first.offset != 0) {
    return 0;
  }
  unsigned shared = seq.shared_prefix_length(first.sequence);
  CHECK_LE(shared, end_pos());
  return shared;
}

scaffold scaffold::rev_comp() const {
  scaffold result = *this;
  result.reverse_in_place();
  return result;
}

void scaffold::reverse_in_place() {
  for (auto& part : m_extents) {
    part.offset = end_pos() - part.offset - part.sequence.size();
    part.sequence = part.sequence.rev_comp();
  }
  std::reverse(m_extents.begin(), m_extents.end());
}

scaffold::iterator scaffold::begin() const {
  if (m_extents.empty()) {
    return end();
  }
  return iterator(this, m_extents.begin(), m_extents.front().sequence.begin(),
                  m_extents.front().offset);
}

scaffold::iterator scaffold::end() const {
  return iterator(this, m_extents.end(), dna_const_iterator(), end_pos());
}

std::pair<dna_slice, dna_slice> scaffold::split_extent_at(aoffset_t start) const {
  for (const auto& part : m_extents) {
    if (part.offset > start) {
      break;
    }
    if (part.offset + aoffset_t(part.sequence.size()) < start) {
      continue;
    }

    aoffset_t part_offset = start - part.offset;
    return std::make_pair(part.sequence.subseq(0, part_offset),
                          part.sequence.subseq(part_offset, part.sequence.size() - part_offset));
  }
  return std::make_pair(dna_slice(), dna_slice());
}

void scaffold::iterator::skip_to(aoffset_t target, const char* description) {
  if (target < m_offset) {
    SPLOG("Scaffold seek error: Seeking to %d from %d, source = %s", target, m_offset, description);
    if (m_scaffold_it == m_scaffold->extents().end()) {
      SPLOG("scaffold iterator at end");
    } else {
      SPLOG("scaffold iterator at scaffold #%ld, pos %d len=%ld",
            m_scaffold_it - m_scaffold->extents().begin(), m_offset,
            m_scaffold_it->sequence.size());
    }
    SPLOG("Scaffolds up to end=%d:", m_scaffold->m_end_pos);
    int n = 0;
    for (const auto& s : m_scaffold->extents()) {
      SPLOG("#%d: %d + %ld -> %ld", n++, s.offset, s.sequence.size(), s.offset + s.sequence.size());
    }
    LOG(FATAL) << "Scaffold seek error; seeking to " << target << " from " << m_offset << " by "
               << description << "\n";
  }

  while (m_offset < target) {
    aoffset_t advance = target - m_offset;

    aoffset_t left_in_this_extent = m_scaffold_it->sequence.end() - m_extent_it;
    if (advance > left_in_this_extent) {
      advance = left_in_this_extent;
    }

    CHECK_GT(advance, 0);

    m_extent_it += advance;
    m_offset += advance;
    if (m_extent_it == m_scaffold_it->sequence.end()) {
      ++m_scaffold_it;
      if (m_scaffold_it == m_scaffold->extents().end()) {
        m_offset = m_scaffold->end_pos();
        return;
      } else {
        m_offset = m_scaffold_it->offset;
        m_extent_it = m_scaffold_it->sequence.begin();
      }
    }
  }
}

}  // namespace variants
