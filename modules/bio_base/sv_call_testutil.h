#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/seq_position.h"
#include "modules/bio_base/sv_call.h"

namespace sv_call_testutil {

// Matches a SV call with the given alleles.  For instance, this
// matches a hetrozygous call on extent #1 at position 1234 facing forward where
// the
// first call is "GATC" and the second is "AATA":
//
// sv_call c;
// EXPECT_THAT(c, SvCall(SeqPosition(1, 1234),
//                       ElementsAre(
//                         Allele("GATC"),
//                         Allele("AATA"))));

class SeqPositionImpl
    : public ::testing::MatcherInterface<const seq_position&> {
 public:
  SeqPositionImpl(int scaffold_id, size_t pos)
      : m_scaffold_id(scaffold_id), m_pos(pos) {}

  bool MatchAndExplain(
      const seq_position& x,
      ::testing::MatchResultListener* listener) const override {
    bool matched = true;
    if (x.scaffold_id != m_scaffold_id) {
      *listener << " where the scaffold id doesn't match";
      matched = false;
    }

    if (x.position != m_pos) {
      *listener << " Where the position doesn't match";
      matched = false;
    }

    return matched;
  };

  void DescribeTo(::std::ostream* os) const override {
    (*os) << "at " << m_pos << " on scaffold " << m_scaffold_id;
  }

 private:
  int m_scaffold_id;
  size_t m_pos;
};

inline ::testing::Matcher<const seq_position&> SeqPosition(int scaffold_id,
                                                           size_t pos) {
  return MakeMatcher(new SeqPositionImpl(scaffold_id, pos));
}

class AlleleImpl : public ::testing::MatcherInterface<const allele&> {
 public:
  AlleleImpl(dna_sequence sequence) : m_sequence(sequence) {}

  bool MatchAndExplain(
      const allele& x,
      ::testing::MatchResultListener* listener) const override {
    return x.seq == m_sequence;
  }

  void DescribeTo(::std::ostream* os) const override {
    (*os) << "has sequence " << m_sequence;
  }

 private:
  dna_sequence m_sequence;
};

inline const ::testing::Matcher<const allele&> Allele(dna_slice slice) {
  return MakeMatcher(new AlleleImpl(dna_sequence(slice.begin(), slice.end())));
}

inline const ::testing::Matcher<const allele&> Allele(const std::string& seq) {
  return MakeMatcher(new AlleleImpl(dna_sequence(seq)));
}

class SvCallImpl : public ::testing::MatcherInterface<const sv_call&> {
 public:
  SvCallImpl(
      const ::testing::Matcher<const seq_position&>& pos_matcher,
      const ::testing::Matcher<const std::vector<allele>&>& alleles_matcher)
      : m_pos_matcher(pos_matcher), m_alleles_matcher(alleles_matcher){};
  bool MatchAndExplain(
      const sv_call& x,
      ::testing::MatchResultListener* listener) const override {
    bool matched = true;
    if (!m_pos_matcher.MatchAndExplain(x.position, listener)) {
      matched = false;
    }

    if (!m_alleles_matcher.MatchAndExplain(x.alleles, listener)) {
      matched = false;
    }
    return matched;
  }

  void DescribeTo(::std::ostream* os) const override {
    (*os) << "is a sv_call ";
    m_pos_matcher.DescribeTo(os);
    (*os) << " containing alleles where ";
    m_alleles_matcher.DescribeTo(os);
  }

 private:
  ::testing::Matcher<const seq_position&> m_pos_matcher;
  ::testing::Matcher<const std::vector<allele>&> m_alleles_matcher;
};

inline const ::testing::Matcher<const sv_call&> SvCall(
    const ::testing::Matcher<const seq_position&>& pos_matcher,
    const ::testing::Matcher<const std::vector<allele>&>& alleles_matcher) {
  return MakeMatcher(new SvCallImpl(pos_matcher, alleles_matcher));
}

std::ostream& operator<<(std::ostream& os, const seq_position& pos) {
  os << "scaffold " << pos.scaffold_id << " position " << pos.position;
  return os;
}

std::ostream& operator<<(std::ostream& os, const allele& a) {
  os << a.seq;
  if (!a.sub_ids.empty()) {
    os << " ids: {";
    for (size_t i = 0; i < a.sub_ids.size(); ++i) {
      if (i != 0) {
        os << ", ";
      }
      os << (a.sub_ids[i]/2) << "*2";
      if (a.sub_ids[i]&1) {
        os << "+1";
      }
    }
    os << "}";
  }
  if (a.depth.size() < 100) {
    os << " depth: {";
    for (size_t i = 0; i < a.depth.size(); ++i) {
      if (i > 0) {
        os << ", ";
      }
      os << a.depth[i];
    }
    os << "}";
  }
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const struct_var& sv) {
  CHECK_GE(sv.var_end, sv.var_start);
  int64_t var_size = sv.var_end - sv.var_start;
  os << "[id=" << sv.var_id << " " << (sv.is_structural ? "SV" : "non-SV")
     << " " << sv.ref_start << (sv.rev_start ? "(RC)" : "") << " " << sv.ref_end
     << (sv.rev_end ? "(RC)" : "") << " var size=" << var_size;
  if (sv.ref_start.scaffold_id == sv.ref_end.scaffold_id) {
    int64_t ref_size = labs(int64_t(sv.ref_end.position) - int64_t(sv.ref_start.position));
    int64_t size_diff = var_size - ref_size;
    os << " ref size: " << ref_size <<" size diff: " << size_diff;
  } else {
    os << " chross-chromosome";
  }
  os << " var: " << sv.assembled.subseq(sv.var_start, sv.var_end - sv.var_start)
     << " assembly: " << sv.assembled.subseq(0, sv.var_start) << " [["
     << sv.assembled.subseq(sv.var_start, sv.var_end - sv.var_start) << "]] "
     << sv.assembled.subseq(sv.var_end, sv.assembled.size() - sv.var_end)
     << " avg depth: " << sv.avg_depth << "]";
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const sv_call& call) {
  os << "sv_call at " << call.position << " with " << call.alleles.size()
     << " alleles:\n";
  for (size_t i = 0; i < call.alleles.size(); ++i) {
    os << " allele " << i << ": " << call.alleles[i] << "\n";
  }
  os << "and " << call.sources.size() << " sources:\n";
  for (size_t i = 0; i < call.sources.size(); ++i) {
    if (call.sources[i].is_structural) {
      os << " source " << i << ": " << call.sources[i] << "\n";
    }
  }
  return os;
}

}  // namespace testing
