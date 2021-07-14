#include "modules/bio_base/dna_testutil.h"

#include <gtest/gtest.h>
#include <iostream>
#include <set>

class dna_testutil_environment : public ::testing::Environment {
  void SetUp() override { enable_test_sequence_expansion(); }
};

::testing::Environment* const dna_testutil_env __attribute__((unused)) =
    ::testing::AddGlobalTestEnvironment(new dna_testutil_environment);

dna_sequence dna_test_sequence(const std::string& str) {
  dna_sequence output;

  for (unsigned char c : str) {
    // Each test sequence is bounded by 'C' if forward, and
    // 'G' if the reverse complement.
    output += dna_base('C');

    for (int i = 0; i < (k_dna_test_sequence_length - 2); ++i) {
      if (c & (1 << i)) {
        output += dna_base('T');
      } else {
        output += dna_base('A');
      }
    }

    output += dna_base('C');
  }
  return output;
}

std::ostream& expand_test_sequence(::std::ostream& os, const dna_slice& seq) {
  os << "\"" << seq.as_string() << "\"";

  enum { NON_DECODING, DECODING_FORWARD, DECODING_RC } state = NON_DECODING;
  std::string all_decoded;
  std::string decoded;

  auto flush = [&]() {
    if (all_decoded.empty()) {
      all_decoded = " (";
    } else {
      all_decoded += " + ";
    }
    switch (state) {
      case NON_DECODING:
        all_decoded += "\"" + decoded + "\"";
        break;
      case DECODING_FORWARD:
        all_decoded += "tseq(\"" + decoded + "\")";
        break;
      case DECODING_RC:
        all_decoded += "tseq_rc(\"" + decoded + "\")";
        break;
    }
    decoded.clear();
  };

  auto get_bits = [](const dna_slice& seq, int start_pos, char& res) -> bool {
    res = 0;
    for (int i = 0; i < 8; ++i) {
      if (seq[start_pos + i] == dna_base('A')) {
        // Decode a 0 bit; do nothing.
      } else if (seq[start_pos + i] == dna_base('T')) {
        res |= (1 << i);
      } else {
        // Not in the format we expect.
        return false;
      }
    }
    return true;
  };

  unsigned i = 0;
  while (i < seq.size()) {
    if ((i + k_dna_test_sequence_length) <= seq.size()) {
      // Try to decode a test sequence
      if (seq[i] == dna_base('C') && seq[i + k_dna_test_sequence_length - 1] == dna_base('C')) {
        char ch;
        if (get_bits(seq, i + 1, ch)) {
          if (i > 0 && state != DECODING_FORWARD) {
            flush();
          }
          state = DECODING_FORWARD;
          decoded.append(std::string({ch}));
          i += k_dna_test_sequence_length;
          continue;
        }
      } else if (seq[i] == dna_base('G') &&
                 seq[i + k_dna_test_sequence_length - 1] == dna_base('G')) {
        dna_slice rc = seq.subseq(i, k_dna_test_sequence_length).rev_comp();
        CHECK(rc[0] == dna_base('C') && rc[k_dna_test_sequence_length - 1] == dna_base('C'));
        char ch;
        if (get_bits(rc, 1, ch)) {
          if (i > 0 && state != DECODING_RC) {
            flush();
          }
          state = DECODING_RC;
          decoded = std::string({ch}) + decoded;
          i += k_dna_test_sequence_length;
          continue;
        }
      }
    }

    // No test sequence matched at current position.
    if (i > 0 && state != NON_DECODING) {
      flush();
    }
    state = NON_DECODING;
    decoded.append(seq.subseq(i, 1).as_string());
    i++;
  }
  if (state != NON_DECODING || !all_decoded.empty()) {
    flush();
    os << all_decoded << ")";
  }
  return os;
}

void PrintTo(const dna_slice& seq, ::std::ostream* os) {
  *os << "\n  ";
  expand_test_sequence(*os, seq);
}

void PrintTo(const dna_sequence& seq, ::std::ostream* os) {
  *os << "\n  ";
  expand_test_sequence(*os, dna_slice(seq.begin(), seq.size()));
}

dna_sequence rand_dna_sequence(std::mt19937& random_source, unsigned seq_len) {
  std::uniform_int_distribution<int> base_gen(0, 3);
  dna_sequence seq;

  for (size_t j = 0; j < seq_len; j++) {
    seq.push_back(dna_base(base_gen(random_source)));
  }
  return seq;
}

// Enable expansion of test sequences when writing dna_sequences to a
// std::ostream.
void enable_test_sequence_expansion() { dna_testutil::g_dna_printer = expand_test_sequence; }

// Disable expansion of test sequences when writing dna_sequences to a
// std::ostream.
void disable_test_sequence_expansion() {
  dna_testutil::g_dna_printer = dna_testutil::default_dna_printer;
}

MATCHER_P2(SequenceMatches, seq_regex, size, "") {
  return arg.size() == size && boost::regex_match(arg.as_string(), seq_regex);
}

testing::Matcher<dna_sequence> dna_sequence_matcher::matcher() const {
  return SequenceMatches(m_regex, m_size);
}

static std::map<dna_sequence, std::string> g_print_seq_annotations;

static void add_print_seq_annotation_internal(const dna_sequence& seq, const std::string& description) {
  g_print_seq_annotations.emplace(seq,description);
}
void add_print_seq_annotation(const dna_sequence& seq, const std::string& description) {
  add_print_seq_annotation_internal(seq, "FWD:" + description);
  add_print_seq_annotation_internal(seq.rev_comp(), "REV:" + description);
}

void add_print_seq_annotation(const std::string& seq_str, const std::string& description) {
  add_print_seq_annotation(dna_sequence(seq_str), description);
}

std::ostream& annotated_dna_printer(std::ostream& os, const dna_slice& seq) {
  std::multimap<size_t /* position */, std::string /* string to output */> ann_closes;
  static const char k_start_red[] = "\x1B[31m";
  static const char k_start_blue[] = "\x1B[34m";
  static const char k_reset[] = "\033[0m";

  for (size_t i = 0;; ++i) {
    auto fit = ann_closes.begin();
    while (fit != ann_closes.end()) {
      CHECK_GE(fit->first, i);
      if (fit->first > i) {
        break;
      }
      os << k_start_red << "</" << fit->second << "> ";
      ann_closes.erase(fit);
      fit = ann_closes.begin();
      if (ann_closes.empty()) {
        os << k_reset;
      } else {
        os << k_start_blue;
      }
    }

    if (i == seq.size()) {
      CHECK(ann_closes.empty());
      break;
    }

    for (const auto& ann : g_print_seq_annotations) {
      if (ann.first.size() + i > seq.size()) {
        continue;
      }
      if (seq.subseq(i, ann.first.size()) != ann.first) {
        continue;
      }

      os << k_start_red << " <" << ann.second << ">" << k_start_blue;
      ann_closes.emplace(i + ann.first.size(), ann.second);
    }

    os << seq[i];
  }
  return os;
}

// Enable annotation of DNA sequences
void enable_annotated_sequences() { dna_testutil::g_dna_printer = annotated_dna_printer; }

namespace dna_testutil {

dna_sequence dna_A("A");
dna_sequence dna_C("C");
dna_sequence dna_G("G");
dna_sequence dna_T("T");

dna_sequence drop_front(unsigned num_bases, const dna_sequence& seq) {
  CHECK_GE(seq.size(), num_bases);
  return seq.subseq(num_bases, seq.size() - num_bases);
}

dna_sequence drop_back(unsigned num_bases, const dna_sequence& seq) {
  CHECK_GE(seq.size(), num_bases);
  return seq.subseq(0, seq.size() - num_bases);
}

dna_sequence_matcher SequenceIs(const std::string& seq_str) {
  return dna_sequence_matcher(seq_str);
}

dna_sequence_matcher LongSequenceMatches(const std::string& seq_regex, size_t size) {
  return dna_sequence_matcher(seq_regex, size);
}

}  // namespace dna_testutil
