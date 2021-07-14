#include "modules/variants/assemble.h"

#include <boost/core/demangle.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <chrono>
#include <ctime>
#include <random>

#include "modules/bio_base/readmap.h"
#include "modules/io/parallel.h"
#include "modules/variants/scaffold.h"

namespace variants {

std::ostream& operator<<(std::ostream& os, const aligned_var& v) {
  return os << "[" << v.left_offset << "+" << (v.right_offset - v.left_offset) << "="
            << v.right_offset << "): " << v.seq;
}

constexpr size_t k_max_output_seq_len = 500;

std::ostream& operator<<(std::ostream& os, const assembly& as) {
  os << "Assembly id=" << as.assembly_id << " " << as.tags;
  if (as.min_overlap) {
    os << " min_overlap=" << as.min_overlap;
  }
  size_t num_ids = 0;
  for (const auto& id : as.merged_assembly_ids) {
    ++num_ids;
    if (num_ids > 5) {
      os << ",...";
      break;
    }
    os << "," << id;
  }
  os << " ";
  as.output_offsets(os);
  as.output_other_info(os);
  os << ": ";
  if (as.seq.size() > k_max_output_seq_len) {
    os << as.seq.subseq(0, k_max_output_seq_len / 2) << "..."
       << as.seq.subseq(as.seq.size() - k_max_output_seq_len / 2, k_max_output_seq_len / 2);
  } else {
    os << as.seq;
  }

  os << " (len=" << as.seq.size();
  if (as.score) {
    os << " score=" << as.score;
  }
  if (as.other_depth) {
    os << " other_depth=" << as.other_depth;
  }
  if (as.other_pair_depth) {
    os << " other_pair_depth=" << as.other_pair_depth;
  }
  if (as.ref_depth) {
    os << " ref_depth=" << as.ref_depth;
  }
  if (as.strand_count) {
    os << " strand_count=" << as.strand_count;
  }
  size_t left_pairs = as.left_pair_matches.size();
  size_t right_pairs = as.right_pair_matches.size();
  size_t tot_pairs = left_pairs + right_pairs;
  if (tot_pairs) {
    os << " pair_matches=" << left_pairs << "+" << right_pairs << "=" << tot_pairs;
  }
  os << ")";
  if (as.edge_coverage) {
    os << " edge_cov(" << *as.edge_coverage << ")";
  }
  return os;
}

void assembly::output_other_info(std::ostream& os) const {
  // if (right_pair_matches || left_pair_matches) {
  //   os << ", right pair count=" << right_pair_matches
  //      << ", left pair count = " << left_pair_matches;
  // }
  if (matches_reference) {
    os << ", (ref)";
  } else {
    os << ", !ref";
  }
}

void assembly::output_offsets(std::ostream& os) const {
  bool both_anchors = right_offset && left_offset;
  os << "[";
  os << left_offset << ":" << left_anchor_len;
  if (both_anchors) {
    os << "+" << (right_offset - left_offset);
  }
  os << "=" << right_offset << ":" << right_anchor_len;
  if (both_anchors) {
    aoffset_t svlen = aoffset_t(seq.size()) - (right_offset - left_offset);
    if (svlen) {
      os << ";svlen=" << svlen;
    }
  }
  os << ")";
}

assemble_options assemble_options::g_defaults;

std::string assemble_stats::as_string() const {
  std::stringstream os;
  os << *this;
  return os.str();
}

assemble_lambda_output::assemble_lambda_output(const std::function<void(assembly_ptr)>& output_f,
                                               const std::string& description)
    : m_output_f(output_f), m_description(description) {}

void assemble_lambda_output::on_assembly(assembly_ptr a) { m_output_f(std::move(a)); }

assemble_lambda_copy::assemble_lambda_copy(const std::function<void(const assembly&)>& copy_f,
                                           pipeline_step_t output, const std::string& description)
    : m_copy_f(copy_f), m_output(std::move(output)), m_description(description) {}

void assemble_lambda_copy::on_assembly(assembly_ptr a) {
  m_copy_f(*a);
  m_output->add(std::move(a));
}

std::pair<assembly_ptr, assembly_ptr> split_assembly(assembly_ptr a, aoffset_t split_pos,
                                                     aoffset_t offset_split_pos) {
  CHECK_GE(offset_split_pos, 0);
  CHECK_LE(offset_split_pos, a->right_offset - a->left_offset);
  aoffset_t left_offset = a->left_offset;
  return split_assembly_absoffset(std::move(a), split_pos, offset_split_pos + left_offset);
}

std::pair<assembly_ptr, assembly_ptr> split_assembly_absoffset(
    assembly_ptr a, aoffset_t split_pos, optional_aoffset abs_offset_split_pos) {
  CHECK_GE(split_pos, 0);
  CHECK_LE(split_pos, a->seq.size());
  bool both_anchors = a->left_offset && a->right_offset;

  assembly_ptr left = make_unique<assembly>(*a);
  assembly_ptr right = make_unique<assembly>(*a);

  left->seqset_entries.clear();
  left->rc_seqset_entries.clear();
  right->seqset_entries.clear();
  right->rc_seqset_entries.clear();

  left->right_offset = right->left_offset = abs_offset_split_pos;

  if (!a->coverage.empty()) {
    left->coverage.resize(split_pos + 1);
    right->coverage.resize(aoffset_t(a->seq.size() + 1) - split_pos);
    for (int i = 0; i < split_pos + 1; ++i) {
      left->coverage[i] = a->coverage[i];
    }
    for (int i = split_pos; i < int(a->seq.size() + 1); ++i) {
      right->coverage[i - split_pos] = a->coverage[i];
    }
  }

  for (auto matches : {&assembly::left_pair_matches, &assembly::right_pair_matches}) {
    (*left).*matches = (*a).*matches;
    (*right).*matches = (*a).*matches;
  }

  left->seq = a->seq.subseq(0, split_pos);
  right->seq = a->seq.subseq(split_pos, aoffset_t(a->seq.size()) - split_pos);

  if (a->left_anchor_len > split_pos) {
    right->left_anchor_len = a->left_anchor_len - split_pos;
    left->left_anchor_len = split_pos;
  } else {
    right->left_anchor_len = 0;
  }

  if (a->right_anchor_len > (aoffset_t(a->seq.size()) - split_pos)) {
    left->right_anchor_len = a->right_anchor_len - (aoffset_t(a->seq.size()) - split_pos);
    right->right_anchor_len = (aoffset_t(a->seq.size()) - split_pos);
  } else {
    left->right_anchor_len = 0;
  }

  CHECK_LE(left->left_anchor_len + left->right_anchor_len, left->seq.size())
      << *a << " split into " << *left << " and  " << *right << " at seq offset " << split_pos
      << " and ref offset " << abs_offset_split_pos;
  if (both_anchors) {
    CHECK_LE(left->left_anchor_len + left->right_anchor_len, left->right_offset - left->left_offset)
        << *a << " split into " << *left << " and " << *right << " at seq offset " << split_pos
        << " and ref offset " << abs_offset_split_pos;
  }

  CHECK_LE(right->left_anchor_len + right->right_anchor_len, right->seq.size())
      << *a << " split into " << *left << " and  " << *right << " at seq offset " << split_pos
      << " and ref offset " << abs_offset_split_pos;
  if (both_anchors) {
    CHECK_LE(right->left_anchor_len + right->right_anchor_len,
             right->right_offset - right->left_offset)
        << *a << " split into " << *left << " and  " << *right << " at seq offset " << split_pos
        << " and ref offset " << abs_offset_split_pos;
  }

  left->aligned_variants.clear();
  right->aligned_variants.clear();
  aoffset_t seq_offset = 0;
  if (!a->aligned_variants.empty()) {
    aoffset_t ref_offset = a->left_offset;
    for (const auto& v : a->aligned_variants) {
      seq_offset += v.left_offset - ref_offset;
      ref_offset = v.left_offset;
      aoffset_t seq_offset_end = seq_offset + v.seq.size();

      if (seq_offset >= aoffset_t(left->seq.size()) && v.left_offset >= right->left_offset) {
        right->aligned_variants.push_back(v);
      } else {
        CHECK_LE(seq_offset_end, left->seq.size())
            << dump_assembly_and_vars(*a) << " split into " << *left << " and  " << *right
            << " at seq offset " << split_pos << " and ref offset " << abs_offset_split_pos;
        CHECK_LE(v.right_offset, left->right_offset)
            << dump_assembly_and_vars(*a) << " split into " << *left << " and  " << *right
            << " at seq offset " << split_pos << " and ref offset " << abs_offset_split_pos;
        left->aligned_variants.push_back(v);
      }

      seq_offset += v.seq.size();
      ref_offset = v.right_offset;
    }
  }

  return std::make_pair(std::move(left), std::move(right));
}

void pad_assembly(assembly* a, aoffset_t new_left_offset, aoffset_t new_right_offset,
                  const assemble_options& options) {
  CHECK(a->coverage.empty()) << "Padding with coverage is unimplemented: " << *a;
  if (new_left_offset < a->left_offset) {
    CHECK(options.scaffold);
    scaffold pad_left_s =
        options.scaffold->subscaffold(new_left_offset, a->left_offset - new_left_offset);
    if (pad_left_s.is_simple()) {
      // TODO(nils): Handle the non-simple case
      dna_slice pad_left = pad_left_s.get_simple();
      a->left_anchor_len += pad_left.size();
      a->left_offset -= pad_left.size();
      a->seq = dna_sequence(pad_left.begin(), pad_left.end()) + a->seq;
    }
  }

  if (new_right_offset > a->right_offset) {
    CHECK(options.scaffold);
    scaffold pad_right_s =
        options.scaffold->subscaffold(a->right_offset, new_right_offset - a->right_offset);
    if (pad_right_s.is_simple()) {
      // TODO(nils): Handle the non-simple case
      dna_slice pad_right = pad_right_s.get_simple();
      a->right_anchor_len += pad_right.size();
      a->right_offset += pad_right.size();
      a->seq += pad_right;
    }
  }
}

bool assemble_pipeline_interface::g_verify_order = false;

void assemble_pipeline_interface::add(assembly_ptr a) {
  if (assembly_needs_trace(*a)) {
    std::cout << "IN:  " << description() << " received " << (void*)a.get() << ": " << *a << "\n";
  }
  if (g_verify_order) {
    if (m_expected_order) {
      if (m_last_assembly) {
        CHECK(!m_expected_order(*a, *m_last_assembly))
            << description() << ": Should not have seen " << *m_last_assembly << " before " << *a;
      }
      m_last_assembly = make_unique<assembly>(*a);
    }
    check_assembly(*a, description());
  }
  on_assembly(std::move(a));
}

std::string assemble_pipeline_interface::description() const {
  std::string cls = boost::core::demangle(typeid(*this).name());
  if (cls.size() > 10 && cls.substr(0, 10) == "variants::") {
    cls = cls.substr(10);
  }
  return cls;
}

sorted_output_pipeline_step::~sorted_output_pipeline_step() {
  flush_sorted();
  CHECK(m_output_queue.empty());
  CHECK(m_left_offsets.empty());
}

void sorted_output_pipeline_step::flush_sorted() {
  flush_sorted_to(std::numeric_limits<aoffset_t>::max());
  m_output->flush();
}

void sorted_output_pipeline_step::track_left_offset(aoffset_t offset) {
  CHECK_GE(offset, m_flush_point);
  m_left_offsets.insert(offset);
}

void sorted_output_pipeline_step::untrack_left_offset(aoffset_t offset) {
  auto it = m_left_offsets.find(offset);
  CHECK(it != m_left_offsets.end());
  m_left_offsets.erase(it);
}

void sorted_output_pipeline_step::flush_sorted_to(aoffset_t flush_offset) {
  CHECK_GE(flush_offset, m_flush_point);
  if (!m_left_offsets.empty()) {
    aoffset_t first_tracked = *m_left_offsets.begin();
    if (flush_offset > first_tracked) {
      flush_offset = first_tracked;
    }
  }

  while (!m_output_queue.empty()) {
    auto it = m_output_queue.begin();
    assembly* aptr = it->first;
    if (min(aptr->left_offset, aptr->right_offset) >= flush_offset) {
      return;
    }
    assembly_ptr a = std::move(it->second);
    m_output_queue.erase(it);
    CHECK_EQ(a.get(), aptr);
    m_output->add(std::move(a));
  }

  m_flush_point = flush_offset;
}

void sorted_output_pipeline_step::sort_and_output(assembly_ptr a) {
  if (assembly_needs_trace(*a)) {
    std::cout << "OUT: " << description() << " produced " << (void*)a.get() << ": " << *a << "\n";
  }
  aoffset_t left_offset = min(a->left_offset, a->right_offset);
  CHECK_GE(left_offset, m_flush_point);
  assembly* aptr = a.get();
  m_output_queue.emplace(aptr, std::move(a));
}

void check_assembly_from_user(const assembly& a) {
  if (a.right_offset && a.left_offset) {
    if (a.right_offset < a.left_offset) {
      throw(io_exception(printstring("Right offset must occur after the left offset in assembly %s",
                                     dump_assembly_and_vars(a).c_str())));
    }
    if (a.seq.size() == 0 && a.left_offset == a.right_offset) {
      throw(io_exception(
          printstring("Assembly must not be empty: %s", dump_assembly_and_vars(a).c_str())));
    }
  }
}

void check_assembly(const assembly& a, std::string where) {
  if (a.matches_reference) {
    CHECK_EQ(a.left_anchor_len, 0) << where << ": " << a;
    CHECK_EQ(a.right_anchor_len, 0) << where << ": " << a;
    CHECK_EQ(a.right_offset - a.left_offset, a.seq.size()) << where << ": " << a;
  }

  if (!a.coverage.empty()) {
    CHECK_EQ(a.coverage.size(), a.seq.size() + 1);
  }

  CHECK(a.left_offset || a.right_offset);
  if (!(a.left_offset && a.right_offset)) {
    // TODO(nils): Maybe do some other kind of checks on half-aligned asemblies too
    CHECK(!a.matches_reference);
    CHECK_GT(a.seq.size(), 0);
    return;
  }

  CHECK_LE(a.left_anchor_len, a.seq.size()) << where << ": " << a;
  CHECK_LE(a.right_anchor_len, a.seq.size()) << where << ": " << a;

  CHECK_GE(a.right_offset, a.left_offset) << where << ": " << a;
  CHECK_LE(a.left_anchor_len, a.right_offset - a.left_offset) << where << ": " << a;
  CHECK_LE(a.right_anchor_len, a.right_offset - a.left_offset) << where << ": " << a;

  if (a.left_anchor_len + a.right_anchor_len > a.right_offset - a.left_offset) {
    // Must match reference exactly.
    CHECK_EQ(a.seq.size(), a.right_offset - a.left_offset) << where << ": " << a;
  }
  CHECK(a.seq.size() > 0 || a.right_offset > a.left_offset) << where << ": " << a;

  if (!a.aligned_variants.empty()) {
    aoffset_t ref_offset = a.left_offset;
    aoffset_t seq_offset = 0;
    for (const auto& var : a.aligned_variants) {
      CHECK_GE(var.left_offset, a.left_offset + a.left_anchor_len)
          << where << ": " << dump_assembly_and_vars(a);
      CHECK_GE(var.left_offset, ref_offset) << where << ": " << dump_assembly_and_vars(a);
      CHECK(var.seq.size() > 0 || var.right_offset > var.left_offset)
          << where << ": " << dump_assembly_and_vars(a);
      CHECK_GE(var.right_offset, var.left_offset) << where << ": " << dump_assembly_and_vars(a);
      seq_offset += (var.left_offset - ref_offset);
      ref_offset = var.left_offset;
      CHECK_LE(seq_offset, a.seq.size()) << where << ": " << dump_assembly_and_vars(a);
      aoffset_t var_len = var.seq.size();
      CHECK_LE(seq_offset + var_len, a.seq.size()) << where << ": " << dump_assembly_and_vars(a);
      CHECK_EQ(a.seq.subseq(seq_offset, var_len), var.seq);
      seq_offset += var_len;
      ref_offset += var.right_offset - var.left_offset;
    }
    CHECK_LE(ref_offset, a.right_offset - a.right_anchor_len)
        << where << ": " << dump_assembly_and_vars(a);
    aoffset_t end_ref_len = a.right_offset - ref_offset;
    seq_offset += end_ref_len;
    ref_offset += end_ref_len;
    CHECK_EQ(seq_offset, a.seq.size()) << where << ": " << dump_assembly_and_vars(a);
    CHECK_EQ(ref_offset, a.right_offset) << where << ": " << dump_assembly_and_vars(a);
  }
}

constexpr bool k_merge_debug = false;

assembly_ptr merge_assemblies(const assembly& a, const assembly& b) {
  CHECK(a.coverage.empty()) << "Merging coverage not supported yet";
  CHECK(b.coverage.empty()) << "Merging coverage not supported yet";

  CHECK(a.matches_reference || !a.aligned_variants.empty()) << "Assembly must be aligned: " << a;
  CHECK(b.matches_reference || !b.aligned_variants.empty()) << "Assembly must be aligned: " << b;

  check_assembly(a, "merge_input1");
  check_assembly(b, "merge_input2");

  bool merge_debug = k_merge_debug;
  bool a_is_first = a.left_offset < b.left_offset;
  const assembly& left = a_is_first ? a : b;
  const assembly& right = a_is_first ? b : a;

  if (merge_debug) {
    std::cout << "Attempting to merge:\n";
    std::cout << "Left: " << dump_assembly_and_vars(left) << "\n";
    std::cout << "Right: " << dump_assembly_and_vars(right) << "\n";
  }

  aoffset_t left_seq_size = left.seq.size();
  aoffset_t right_seq_size = right.seq.size();

  if (left.right_offset < right.left_offset) {
    // Disjoint; can't join
    return nullptr;
  }

  assembly_ptr result = make_unique<assembly>();
  *result = left;

  if (right.ml_features || left.ml_features) {
    CHECK(right.ml_features)
        << "Either both assemblies or neither assembly should have ML features present";
    CHECK(left.ml_features)
        << "Either both assemblies or neither assembly should have ML features present";
    if (right.ml_features->alt_seq.size() > left.ml_features->alt_seq.size()) {
      result->ml_features = right.ml_features;
    } else {
      result->ml_features = left.ml_features;
    }
  }

  if (left.right_offset > right.right_offset) {
    // Completely subsumed.
    if (merge_debug) {
      std::cout << "Subsumed\n";
    }
  } else {
    result->right_offset = right.right_offset;
  }

  result->left_anchor_len = left.left_anchor_len;
  if (left.right_offset < right.right_offset) {
    result->right_anchor_len = right.right_anchor_len;
  }

  result->merged_assembly_ids.push_back(right.assembly_id);
  for (const auto& assembly_id : right.merged_assembly_ids) {
    result->merged_assembly_ids.push_back(assembly_id);
  }

  for (const auto& read_id : right.rc_read_ids) {
    result->rc_read_ids.insert(read_id);
  }

  result->score = std::min(left.score, right.score);
  result->matches_reference = left.matches_reference && right.matches_reference;

  for (auto matches : {&assembly::left_pair_matches, &assembly::right_pair_matches}) {
    auto& pm = (*result).*matches;
    const auto& right_pm = right.*matches;
    pm.insert(pm.end(), right_pm.begin(), right_pm.end());
    std::sort(pm.begin(), pm.end());
    pm.erase(std::unique(pm.begin(), pm.end()), pm.end());
  }

  auto left_vit = left.aligned_variants.begin();
  auto right_vit = right.aligned_variants.begin();

  aoffset_t ref_pos = left.left_offset;
  aoffset_t left_seq_pos = 0;
  aoffset_t right_seq_pos = 0;

  auto advance = [&](aoffset_t tot_ref_advance, aoffset_t tot_seq_advance) -> bool {
    if (merge_debug) {
      std::cout << "Advancing ref+" << tot_ref_advance << " seq+" << tot_seq_advance << "\n";
    }
    CHECK_GE(tot_ref_advance, 0);

    aoffset_t ref_advance_remaining = tot_ref_advance;
    aoffset_t seq_advance_remaining = tot_seq_advance;

    while (ref_advance_remaining || seq_advance_remaining) {
      aoffset_t ref_adv = ref_advance_remaining;

      // If we're transitioning between assemblies, only advance part at once.
      if (ref_pos < right.left_offset && ref_pos + ref_adv > right.left_offset) {
        // Only advance until the start of the right sequence.
        ref_adv = right.left_offset - ref_pos;
      }

      if (ref_pos < left.right_offset && ref_pos + ref_adv > left.right_offset) {
        // Only advance until the end of the left sequence
        ref_adv = left.right_offset - ref_pos;
      }

      if (ref_pos < right.right_offset && ref_pos + ref_adv > right.right_offset) {
        // Only advance until the end of the right sequence (right
        // sequence is subsumed in left sequence)
        ref_adv = right.right_offset - ref_pos;
      }

      bool left_finished = ref_pos >= left.right_offset;
      bool right_started = ref_pos >= right.left_offset;
      bool right_finished = ref_pos >= right.right_offset;

      aoffset_t seq_adv = seq_advance_remaining;

      if (ref_pos == left.right_offset && left_seq_pos < left_seq_size) {
        // Insert on left side
        left_finished = false;
        if (ref_pos == right.left_offset) {
          right_started = false;
        }
        seq_adv = std::min(seq_adv, left_seq_size - left_seq_pos);
        ref_adv = 0;

        if (!seq_adv && !ref_adv) {
          if (merge_debug) {
            std::cout << "No seq advance, but still seq left to advance through\n";
          }
          return false;
        }
      }

      if (ref_adv != ref_advance_remaining || seq_adv != seq_advance_remaining) {
        if (tot_ref_advance != tot_seq_advance) {
          if (merge_debug) {
            std::cout << "Spans assembly boundaries, but size mismatch\n";
          }
          return false;
        }
        seq_adv = ref_adv;
        if (!seq_adv && !ref_adv) {
          if (merge_debug) {
            std::cout << "No seq advance, but still seq left to advance through(2)\n";
          }
          return false;
        }
      }

      if (left_finished) {
        if (left_seq_pos != left_seq_size) {
          if (merge_debug) {
            std::cout << "Not all of left consumed.\n";
          }
          return false;
        }
      }
      if (right_finished) {
        if (right_seq_pos != right_seq_size) {
          if (merge_debug) {
            std::cout << "Not all of right consumed.\n";
          }
          return false;
        }
        CHECK(right_started);
      }

      CHECK_GE(ref_adv, 0);

      if (merge_debug) {
        std::cout << "Advancing partially, " << ref_pos << "+" << ref_adv
                  << ", seq_adv = " << seq_adv << "\n";
        std::cout << "Left remaining:  "
                  << left.seq.subseq(left_seq_pos, left_seq_size - left_seq_pos) << "\n";
        std::cout << "Right remaining: "
                  << right.seq.subseq(right_seq_pos, right_seq_size - right_seq_pos) << "\n";
      }

      CHECK(seq_adv || ref_adv) << "seq adv: " << seq_adv << " ref adv: " << ref_adv;

      seq_advance_remaining -= seq_adv;
      ref_advance_remaining -= ref_adv;
      CHECK_GE(seq_advance_remaining, 0);
      CHECK_GE(ref_advance_remaining, 0);

      if (!left_finished && !right_started) {
        CHECK(!right_finished);
        if (merge_debug) {
          std::cout << "Left unique: " << left.seq.subseq(left_seq_pos, seq_adv) << "\n";
        }
        // In the portion of the left that doesn't overlap the right.  Only advance the left.
        ref_pos += ref_adv;
        if (left_seq_pos + seq_adv > left_seq_size) {
          if (merge_debug) {
            std::cout << "Ran out of seq in left unique region\n";
          }
          return false;
        }
        left_seq_pos += seq_adv;

        CHECK_LE(ref_pos, left.right_offset);
        CHECK_LE(ref_pos, right.left_offset);
        CHECK_LE(left_seq_pos, left_seq_size) << "a: " << a << " b: " << b;
        CHECK_EQ(right_seq_pos, 0);
      } else if (right_started && !left_finished && !right_finished) {
        // Overlapping region where we're in both left and right
        CHECK_LE(ref_pos + ref_adv, left.right_offset) << "a: " << a << " b: " << b;
        CHECK_LE(ref_pos + ref_adv, right.right_offset);
        if (left_seq_pos + seq_adv > left_seq_size || right_seq_pos + seq_adv > right_seq_size) {
          if (merge_debug) {
            std::cout << "Ran out of seq in shared region\n";
          }
          return false;
        }
        dna_slice left_shared = dna_slice(left.seq).subseq(left_seq_pos, seq_adv);
        dna_slice right_shared = dna_slice(right.seq).subseq(right_seq_pos, seq_adv);
        if (left_shared != right_shared) {
          if (merge_debug) {
            std::cout << "left shared " << left_shared << " != right shared " << right_shared
                      << "\n";
          }
          return false;
        }
        if (merge_debug) {
          std::cout << "Shared: " << left_shared << "\n";
        }

        ref_pos += ref_adv;
        left_seq_pos += seq_adv;
        right_seq_pos += seq_adv;

        CHECK_LE(ref_pos, right.right_offset);
        CHECK_LE(ref_pos, left.right_offset);
        CHECK_LE(left_seq_pos, left_seq_size);
        CHECK_LE(right_seq_pos, right_seq_size);
      } else if (left_finished) {
        CHECK(right_started);
        CHECK(!right_finished);
        // We're past the overlapping region between left and right, and are in right only.
        if (merge_debug) {
          std::cout << "Right unique: " << right.seq.subseq(right_seq_pos, seq_adv) << "\n";
        }

        result->seq += right.seq.subseq(right_seq_pos, seq_adv);
        ref_pos += ref_adv;
        if (right_seq_pos + seq_adv > right_seq_size) {
          if (merge_debug) {
            std::cout << "Ran out of seq in right unique region\n";
          }
          return false;
        }
        right_seq_pos += seq_adv;

        CHECK_GE(ref_pos, left.right_offset);
        CHECK_LE(ref_pos, right.right_offset);
        CHECK_EQ(left_seq_pos, left_seq_size);
        CHECK_LE(right_seq_pos, right_seq_size);
      } else if (right_finished && !left_finished) {
        CHECK(right_started);
        // Right is completely subsumed by left, and we're past the subsumed region.
        CHECK_GE(ref_pos, right.right_offset);
        CHECK_LT(ref_pos, left.right_offset);

        if (merge_debug) {
          std::cout << "Left unique, after right: " << left.seq.subseq(left_seq_pos, seq_adv)
                    << "\n";
        }

        ref_pos += ref_adv;
        left_seq_pos += seq_adv;

        CHECK_EQ(right_seq_pos, right_seq_size);
        CHECK_LE(left_seq_pos, left_seq_size);
        CHECK_LE(ref_pos, left.right_offset);
        CHECK_GE(ref_pos, right.right_offset);
      } else {
        LOG(FATAL) << "Unable to determine where " << ref_pos << " lies merging " << left << " and "
                   << right << "\n";
      }
      if (merge_debug) {
        std::cout << "Advance complete, ref=" << ref_pos << " left seq=" << left_seq_pos << "/"
                  << left_seq_size << " right seq=" << right_seq_pos << "/" << right_seq_size
                  << "\n";
      }
    }
    return true;
  };

  while (left_vit != left.aligned_variants.end() && right_vit != right.aligned_variants.end()) {
    if (merge_debug) {
      std::cout << "\nComparing aligned, ref_pos=" << ref_pos << ": " << *left_vit << " vs "
                << *right_vit << "\n";
    }
    // Advance reference until the left edge of this variant.
    aoffset_t ref_adv = left_vit->left_offset - ref_pos;
    if (!advance(ref_adv, ref_adv)) {
      if (merge_debug) {
        std::cout << "Ref mismatch\n";
      }
      return nullptr;
    }

    if (left_vit->right_offset > right.left_offset) {
      // Variant should be shared between left and right.
      if (*right_vit != *left_vit) {
        if (merge_debug) {
          std::cout << "Variants don't match: " << *left_vit << " != " << *right_vit << "\n";
        }
        return nullptr;
      }

      if (!advance(left_vit->right_offset - left_vit->left_offset, left_vit->seq.size())) {
        return nullptr;
      }
      CHECK_EQ(ref_pos, left_vit->right_offset);

      ++left_vit;
      ++right_vit;
    } else {
      // Variant in left before right
      if (merge_debug) {
        std::cout << "Variant in left before right; advancing left only\n";
      }

      if (!advance(left_vit->right_offset - left_vit->left_offset, left_vit->seq.size())) {
        return nullptr;
      }
      CHECK_EQ(ref_pos, left_vit->right_offset);
      ++left_vit;
    }
  }

  while (left_vit != left.aligned_variants.end()) {
    // Variant in left after right.
    if (merge_debug) {
      std::cout << "\nUnique left variant, ref_pos=" << ref_pos << ": " << *left_vit << "\n";
    }

    // Advance through reference
    aoffset_t ref_adv = left_vit->left_offset - ref_pos;
    if (!advance(ref_adv, ref_adv)) {
      return nullptr;
    }

    // Advance through variant
    if (!advance(left_vit->right_offset - left_vit->left_offset, left_vit->seq.size())) {
      return nullptr;
    }
    CHECK_EQ(ref_pos, left_vit->right_offset);
    ++left_vit;
  }

  while (right_vit != right.aligned_variants.end()) {
    // Variant in right after left.
    if (merge_debug) {
      std::cout << "\nUnique right variant, ref_pos=" << ref_pos << ": " << *right_vit << "\n";
    }

    // Advance through reference
    aoffset_t ref_adv = right_vit->left_offset - ref_pos;
    if (!advance(ref_adv, ref_adv)) {
      return nullptr;
    }

    // Advance through variant
    if (!advance(right_vit->right_offset - right_vit->left_offset, right_vit->seq.size())) {
      return nullptr;
    }
    result->aligned_variants.push_back(*right_vit);
    CHECK_EQ(ref_pos, right_vit->right_offset);
    ++right_vit;
  }

  // Advance through reference at end
  aoffset_t last_ref_adv = result->right_offset - ref_pos;
  if (merge_debug) {
    std::cout << "Variants merged; still need to advance " << last_ref_adv << " to "
              << result->right_offset << "\n";
    std::cout << "Left remaining:  " << left.seq.subseq(left_seq_pos, left_seq_size - left_seq_pos)
              << "\n";
    std::cout << "Right remaining: "
              << right.seq.subseq(right_seq_pos, right_seq_size - right_seq_pos) << "\n";
  }
  if (!advance(last_ref_adv, last_ref_adv)) {
    return nullptr;
  }

  if (left_seq_pos != left_seq_size) {
    if (merge_debug) {
      std::cout << "Did not consume all of left\n";
    }
    return nullptr;
  }

  if (right_seq_pos != right_seq_size) {
    if (merge_debug) {
      std::cout << "Did not consume all of right\n";
    }
    return nullptr;
  }

  CHECK_EQ(ref_pos, result->right_offset);

  if (merge_debug) {
    std::cout << "Merge result: " << dump_assembly_and_vars(*result) << "\n";
  }
  check_assembly(*result, "merged_assembly");
  return result;
}

std::string dump_coverage(const std::vector<int>& cov) {
  std::stringstream os;

  for (unsigned idx = 0; idx < cov.size(); ++idx) {
    int c = cov[idx];
    constexpr char covstr[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (c == k_max_output_seq_len && cov.size() > 2 * k_max_output_seq_len) {
      os << "  ...";
    } else if (idx < k_max_output_seq_len || idx + k_max_output_seq_len > cov.size()) {
      if (c >= int(strlen(covstr)) || c < 0) {
        os << " " << c;
      } else {
        os << covstr[c];
      }
    }
  }
  return os.str();
}

std::string dump_assembly_and_vars(const assembly& a) {
  std::stringstream os;
  os << a;
  os << " vars:\n";
  for (const auto& v : a.aligned_variants) {
    os << "  " << v << "\n";
  }
  if (!a.coverage.empty()) {
    os << "Cov: " << dump_coverage(a.coverage) << " (min=" << container_min(a.coverage) << ")\n";
  }
  if (!a.pair_coverage.empty()) {
    os << "Pair:" << dump_coverage(a.pair_coverage) << " (min=" << container_min(a.coverage)
       << ")\n";
  }
  return os.str();
}

size_t allocate_assembly_id() {
  static std::atomic<uint64_t> g_assembly_id{1};

  return g_assembly_id.fetch_add(1);
}

namespace {

std::string print_read_ids(const read_id_set& read_ids) { return std::to_string(read_ids.size()); }

}  // namespace

std::ostream& operator<<(std::ostream& os, const edge_coverage_t& cov) {
  return os << "var_start=" << print_read_ids(cov.variant_start)
            << " var_end=" << print_read_ids(cov.variant_end)
            << " ref_start=" << print_read_ids(cov.reference_start)
            << " ref_end=" << print_read_ids(cov.reference_end)
            << " interior=" << print_read_ids(cov.interior);
}

std::set<size_t> g_trace_assembly_ids;
std::set<aoffset_t> g_trace_offsets;
bool g_trace_all_assemblies = false;

#if NDEBUG
static constexpr bool k_allow_assembly_trace = false;
#else
static constexpr bool k_allow_assembly_trace = true;
#endif

void add_assembly_trace(size_t assembly_id) {
  if (k_allow_assembly_trace) {
    g_trace_assembly_ids.insert(assembly_id);
  } else {
    SPLOG(
        "WARNING: Assembly tracing of assembly_id %ld is incompatible with NDEBUG compilation mode",
        assembly_id);
  }
}

void add_offset_trace(aoffset_t offset) {
  if (k_allow_assembly_trace) {
    g_trace_offsets.insert(offset);
  } else {
    SPLOG("WARNING: Offset tracing of offset %d is incompatible with NDEBUG compilation mode",
          offset);
  }
}

void reset_assembly_trace() {
  g_trace_assembly_ids.clear();
  g_trace_offsets.clear();
}

void reverse_assembly_in_place(assembly* a, const readmap* rm, aoffset_t ref_end_pos) {
  a->seq = a->seq.rev_comp();
  std::swap(a->left_anchor_len, a->right_anchor_len);
  std::swap(a->left_offset, a->right_offset);
  if (a->left_offset) {
    a->left_offset = ref_end_pos - a->left_offset;
  }
  if (a->right_offset) {
    a->right_offset = ref_end_pos - a->right_offset;
  }
  for (auto& v : a->aligned_variants) {
    v.seq = v.seq.rev_comp();
    v.left_offset = ref_end_pos - v.left_offset;
    v.right_offset = ref_end_pos - v.right_offset;
    std::swap(v.left_offset, v.right_offset);
  }
  std::reverse(a->aligned_variants.begin(), a->aligned_variants.end());
  std::swap(a->left_pair_matches, a->right_pair_matches);
  a->seqset_entries.swap(a->rc_seqset_entries);
  if (rm) {
    for (auto* matches : {&a->left_pair_matches, &a->right_pair_matches}) {
      std::vector<uint32_t> old_matches = std::move(*matches);
      matches->clear();
      matches->reserve(old_matches.size());
      for (uint32_t read_id : old_matches) {
        matches->push_back(rm->get_rev_comp(read_id));
      }
    }
    for (auto* cov : {&a->read_coverage, &a->pair_read_coverage}) {
      if (!*cov) {
        continue;
      }

      read_coverage_set new_entries;
      for (const auto& cov_entry : (*cov)->reads()) {
        read_coverage_read_t new_entry;
        new_entry.offset = aoffset_t(a->seq.size()) - cov_entry.read_len - cov_entry.offset;
        new_entry.read_len = cov_entry.read_len;
        for (uint32_t read_id : cov_entry.read_ids) {
          readmap::read rd = rm->get_read_by_id(read_id);
          new_entry.read_ids.insert(rd.get_rev_comp().get_read_id());
        }
        new_entries.insert(std::move(new_entry));
      }
      cov->emplace(new_entries.build_and_clear(a->seq.size()));
    }
    read_id_set new_rc_read_ids;
    for (uint32_t read_id : a->rc_read_ids) {
      new_rc_read_ids.insert(rm->get_rev_comp(read_id));
    }
    a->rc_read_ids = std::move(new_rc_read_ids);

    if (a->edge_coverage) {
      auto& ec = *a->edge_coverage;
      std::swap(ec.variant_start, ec.variant_end);
      std::swap(ec.reference_start, ec.reference_end);

      for (auto* collection : {&ec.variant_start, &ec.variant_end, &ec.interior,
                               &ec.reference_start, &ec.reference_end}) {
        read_id_set reversed;
        for (uint32_t read_id : *collection) {
          reversed.insert(rm->get_rev_comp(read_id));
        }
        std::swap(*collection, reversed);
      }
    }
  }
}

half_aligned_assembly reverse_half_aligned(half_aligned_assembly ha, const readmap* rm,
                                           aoffset_t ref_end_pos) {
  read_id_set reversed_read_ids;
  for (uint32_t read_id : ha.rc_read_ids) {
    reversed_read_ids.insert(rm->get_rev_comp(read_id));
  }
  ha.rc_read_ids = std::move(reversed_read_ids);
  ha.offset = ref_end_pos - ha.offset;
  ha.right_anchor = !ha.right_anchor;
  ha.seq = ha.seq.rev_comp();
  return ha;
}

std::ostream& operator<<(std::ostream& os, const half_aligned_assembly& ha) {
  if (ha.right_anchor) {
    os << "[?, " << ha.scaffold_name << ":" << ha.offset << ")";
  } else {
    os << "[" << ha.scaffold_name << ":" << ha.offset << ", ?)";
  }
  os << " " << ha.seq << " id=" << ha.assembly_id << " reads=" << ha.rc_read_ids.size();
  return os;
}

// Default to GRAPH_DISCOVER unless we need the compatibility sort
// order in 'biograph variants'.
sort_order canon_assembly_order::g_default_sort_order = sort_order::GRAPH_DISCOVER;

bool canon_assembly_order::operator()(const assembly& a, const assembly& b) const {
#define COMPARE(FIELD)        \
  if (a.FIELD != b.FIELD) {   \
    return a.FIELD < b.FIELD; \
  }
#define COMPARE_R(FIELD)      \
  if (a.FIELD != b.FIELD) {   \
    return a.FIELD > b.FIELD; \
  }

  if (m_sort_order == sort_order::LEFT_OFFSET_ONLY) {
    // Old sort order only sorts by left offset.  TODO(nils): Figure out
    // why we need to support old_sort_order=true and make it not
    // required.
    COMPARE(left_offset);
    return false;
  }
  aoffset_t a_min = min(a.left_offset, a.right_offset);
  aoffset_t b_min = min(b.left_offset, b.right_offset);
  if (a_min != b_min) {
    return a_min < b_min;
  }

  COMPARE(matches_reference);

  // fully aligned before half aligned.
  if (bool(a.left_offset) != bool(b.left_offset)) {
    return bool(a.left_offset);
  }
  if (bool(a.right_offset) != bool(b.right_offset)) {
    return bool(a.right_offset);
  }

  aoffset_t a_max = max(a.left_offset, a.right_offset);
  aoffset_t b_max = max(b.left_offset, b.right_offset);
  if (a_max != b_max) {
    return a_max > b_max;
  }

  if (m_sort_order == sort_order::GRAPH_DISCOVER) {
    // Make sure things witth the same location and sequences get sorted
    // together so they can be deduplicated.
    COMPARE(seq);
  }

  COMPARE(tags);
  COMPARE(left_anchor_len);
  COMPARE(right_anchor_len);
  COMPARE(score);
  COMPARE(left_pair_matches.size());
  COMPARE(right_pair_matches.size());
  COMPARE(rc_read_ids.size());
  COMPARE_R(seq);
  return false;
#undef COMPARE
#undef COMPARE_R
}

std::string sorted_output_pipeline_step::sorted_output_stats(
    boost::optional<aoffset_t> relative_to) const {
  std::stringstream out;

  out << "Sorted output flush point = " << m_flush_point;
  if (relative_to) {
    out << "(" << (*relative_to - m_flush_point) << " behind)";
  }
  out << ", " << m_output_queue.size() << " sorted assemblies queued, " << m_left_offsets.size()
      << " left offsets tracked";
  if (!m_left_offsets.empty()) {
    aoffset_t earliest_tracked = *m_left_offsets.begin();
    out << ", earliest tracked=" << earliest_tracked;
    if (relative_to) {
      out << "(" << (*relative_to - earliest_tracked) << " behind)";
    }
  }
  return out.str();
}

std::string string_set::to_string() const {
  bool first = true;
  std::stringstream out;
  out << "(";
  for (const auto& pid : *this) {
    if (!first) {
      out << ",";
    }
    first = false;
    out << pid;
  }
  out << ")";
  return out.str();
}

std::string string_set::to_string_short() const {
  bool first = true;
  std::stringstream out;
  for (const auto& pid : *this) {
    if (!first) {
      out << ",";
    }
    first = false;
    out << pid;
  }
  return out.str();
}

string_set& string_set::operator+=(const string_set& rhs) {
  insert(rhs.begin(), rhs.end());
  return *this;
}

string_set string_set::operator+(const string_set& rhs) const {
  string_set result = *this;
  result += rhs;
  return result;
}

string_set& string_set::operator-=(const string_set& rhs) {
  auto it = begin();
  auto rhs_it = rhs.begin();

  while (it != end() && rhs_it != rhs.end()) {
    if (*it < *rhs_it) {
      ++it;
      continue;
    }
    if (*it > *rhs_it) {
      ++rhs_it;
      continue;
    }
    it = erase(it);
    ++rhs_it;
  }
  return *this;
}

string_set string_set::operator-(const string_set& rhs) const {
  string_set result = *this;
  result -= rhs;
  return result;
}

string_set& string_set::operator&=(const string_set& rhs) {
  auto it = begin();
  auto rhs_it = rhs.begin();

  while (it != end() && rhs_it != rhs.end()) {
    if (*it < *rhs_it) {
      it = erase(it);
      continue;
    }
    if (*it > *rhs_it) {
      ++rhs_it;
      continue;
    }
    ++it;
    ++rhs_it;
  }
  erase(it, end());
  return *this;
}

string_set string_set::operator&(const string_set& rhs) const {
  string_set result;
  auto it = begin();
  auto rhs_it = rhs.begin();

  while (it != end() && rhs_it != rhs.end()) {
    if (*it < *rhs_it) {
      ++it;
      continue;
    }
    if (*it > *rhs_it) {
      ++rhs_it;
      continue;
    }

    result.insert(result.end(), *it);
    ++it;
    ++rhs_it;
  }
  return result;
}

static constexpr size_t k_count_assemblies = 0;
// static constexpr size_t k_count_assemblies = 500000;
static std::atomic<bool> g_sample_assembly{false};

std::atomic<size_t> g_tot_assemblies;
std::atomic<size_t> g_cum_assemblies;

namespace {

void inc_assembly_count() {
  if (!k_count_assemblies) {
    return;
  }
  ++g_tot_assemblies;
  if (((++g_cum_assemblies) % (k_count_assemblies ? k_count_assemblies : 1) == 0)) {
    std::cerr << "asm_alloc(" << g_tot_assemblies << "/" << g_cum_assemblies << "-"
              << g_cum_assemblies - g_tot_assemblies << ")\n";
    g_sample_assembly = true;
  }
}

void dec_assembly_count(assembly* a) {
  if (!k_count_assemblies) {
    return;
  }
  --g_tot_assemblies;
  bool expected = true;
  if (g_sample_assembly.compare_exchange_weak(expected, false)) {
    std::cerr << "Sampled assembly: " << *a << "\n";
  }
}

}  // namespace

assembly::assembly(optional_aoffset left_off, optional_aoffset right_off, dna_sequence aseq)
    : assembly_id(allocate_assembly_id()),
      left_offset(left_off),
      right_offset(right_off),
      seq(aseq) {
  inc_assembly_count();
}

assembly::assembly(optional_aoffset left_off, optional_aoffset right_off, dna_sequence aseq,
                   size_t asm_id)
    : assembly_id(asm_id), left_offset(left_off), right_offset(right_off), seq(aseq) {
  inc_assembly_count();
}

assembly::assembly() { inc_assembly_count(); }

assembly::assembly(const assembly& rhs) {
  inc_assembly_count();
  (*this) = rhs;
}
assembly::assembly(assembly&& rhs) {
  inc_assembly_count();
  (*this) = std::move(rhs);
}
assembly& assembly::operator=(const assembly&) = default;
assembly& assembly::operator=(assembly&&) = default;
assembly::~assembly() { dec_assembly_count(this); }

bool assembly_needs_trace(const assembly& a) {
#if NDEBUG
  constexpr bool k_enable_assembly_tracing = false;
#else
  constexpr bool k_enable_assembly_tracing = true;
#endif
  if (!k_enable_assembly_tracing) {
    return false;
  }

  if (g_trace_all_assemblies) {
    return true;
  }

  if (g_trace_assembly_ids.empty()) {
    return false;
  }

  if (g_trace_assembly_ids.count(a.assembly_id)) {
    return true;
  }
  for (size_t aid : a.merged_assembly_ids) {
    if (g_trace_assembly_ids.count(aid)) {
      return true;
    }
  }
  return false;
}

bool offset_needs_trace(aoffset_t offset) {
#if NDEBUG
  constexpr bool k_enable_assembly_tracing = false;
#else
  constexpr bool k_enable_assembly_tracing = true;
#endif

  if (!k_enable_assembly_tracing) {
    return false;
  }

  if (g_trace_offsets.count(offset)) {
    return true;
  }
  return false;
}

optional_aoffset min(const optional_aoffset& lhs, const optional_aoffset& rhs) {
  if (lhs && rhs) {
    return std::min<aoffset_t>(lhs, rhs);
  } else if (lhs) {
    return lhs;
  } else {
    return rhs;
  }
}

optional_aoffset max(const optional_aoffset& lhs, const optional_aoffset& rhs) {
  if (lhs && rhs) {
    return std::max<aoffset_t>(lhs, rhs);
  } else if (lhs) {
    return lhs;
  } else {
    return rhs;
  }
}

const optional_aoffset optional_aoffset::none;

const absl::btree_set<seqset_range> seqset_path::g_empty;

void seqset_path::propagate_from_end(const absl::btree_set<seqset_range>& new_ends, dna_slice seq,
                                     const assemble_options& opts) {
  constexpr bool k_dbg = false;
  extern std::ostream& operator<<(std::ostream&, const absl::btree_set<seqset_range>&);

  if (k_dbg) {
    std::cerr << "Starting prop from end: " << new_ends << " seq=" << seq << "\n";
  }
  if (!empty() && new_ends == ends()) {
    if (k_dbg) {
      std::cerr << "Skipping, already done\n";
    }
    return;
  }

  absl::btree_set<seqset_range> cur = new_ends;
  seqset_set_dedup_prefixes(cur);
  m_entries[seq.size()] = cur;
  m_entries.try_emplace(0);  // Just in case it doesn't exist already.
  CHECK_EQ(m_entries.begin()->first, 0);
  CHECK_EQ(m_entries.rbegin()->first, seq.size());

  aoffset_t offset = seq.size();
  auto it = m_entries.rbegin();
  CHECK_EQ(it->first, offset);
  CHECK(it != m_entries.rend());
  ++it;

  while (offset) {
    --offset;
    absl::btree_set<seqset_range> new_cur;
    for (const auto& r : cur) {
      new_cur.insert(r.push_front_drop(seq[offset]));
    }
    seqset_set_dedup_prefixes(new_cur);
    cur = std::move(new_cur);

    if (k_dbg) {
      std::cerr << "offset=" << offset << " cur=" << cur << "\n";
    }

    CHECK(it != m_entries.rend());
    if (it->first == offset) {
      it->second = cur;
      ++it;
    } else {
      CHECK_LT(it->first, offset);
    }

    if (opts.readmap) {
      for (const auto& r : cur) {
        unsigned nreads = 0;
        std::vector<seqset_range> found_mates;
        for (const auto& rd : opts.readmap->get_prefix_reads(r)) {
          if (!rd.has_mate()) {
            continue;
          }

          ++nreads;
          if (nreads > opts.max_pairs_per_read) {
            break;
          }
          found_mates.push_back(rd.get_mate_rc().get_seqset_entry());
        }
        if (nreads > opts.max_pairs_per_read) {
          continue;
        }

        m_mates.insert(found_mates.begin(), found_mates.end());
      }
    }
  }
  CHECK(it == m_entries.rend());
}

void seqset_path::clear() {
  m_entries.clear();
  m_mates.clear();
}

const absl::btree_set<seqset_range>& seqset_path::starts() const {
  if (empty()) {
    return g_empty;
  } else {
    return m_entries.begin()->second;
  }
}

const absl::btree_set<seqset_range>& seqset_path::ends() const {
  if (empty()) {
    return g_empty;
  } else {
    return m_entries.rbegin()->second;
  }
}

size_t seqset_path::size() const {
  CHECK(!empty());
  return m_entries.rbegin()->first;
}

absl::btree_set<seqset_range> seqset_path::entries_at_offset(aoffset_t offset) {
  absl::btree_set<seqset_range> entries;
  auto range = m_entries.equal_range(offset);
  for (auto it = range.first; it != range.second; ++it) {
    entries.insert(it->second.begin(), it->second.end());
  }
  return entries;
}

void seqset_path::add(aoffset_t offset, seqset_range r) {
  auto& e = m_entries[offset];
  e.insert(r);
  seqset_set_dedup_prefixes(e);
}
void seqset_path::add(aoffset_t offset, const absl::btree_set<seqset_range>& rs) {
  auto& e = m_entries[offset];
  e.insert(rs.begin(), rs.end());
  seqset_set_dedup_prefixes(e);
}

seqset_path::~seqset_path() {}

void seqset_path::swap(seqset_path& rhs) {
  m_entries.swap(rhs.m_entries);
  m_mates.swap(rhs.m_mates);
}

void seqset_set_dedup_prefixes(absl::btree_set<seqset_range>& rs) {
  if (rs.empty()) {
    return;
  }

  bool again = true;
  while (again) {
    again = false;
    for (auto it = rs.begin(); it != rs.end(); ++it) {
      auto next = it;
      ++next;
      for (; next != rs.end(); ++next) {
        if (next->begin() >= it->end()) {
          break;
        }

        if (it->begin() <= next->begin() && it->end() >= next->end()) {
          // "it" is prefix of "next"; discard it.
          it = rs.erase(it);
          again = true;
          break;
        }
      }
    }
  }
}

}  // namespace variants
