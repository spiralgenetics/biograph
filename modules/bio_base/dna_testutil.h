#pragma once

#include "modules/bio_base/dna_sequence.h"

#include <random>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <boost/regex.hpp>

// Generates unique DNA sequences that can be chained together to
// simulate structures that are uniquely identifiable.
//
// Each character in the input string is translated to a sequence of
// k_dna_test_sequence_length bases that is unique for that character,
// even if its reverse complement is taken.
//
// This allows us to generate simulated reads and explicitly control
// where they match each other.
dna_sequence dna_test_sequence(const std::string& str);

// Returns a random dna sequence of the given length.
dna_sequence rand_dna_sequence(std::mt19937& random_source, unsigned seq_len);

// The length of a sequence returned by dna_test_sequence.
constexpr int k_dna_test_sequence_length = 10;

// Allow gtest to pretty print values of type 'dna_sequence':
void PrintTo(const dna_slice& seq, ::std::ostream* os);
void PrintTo(const dna_sequence& seq, ::std::ostream* os);

// Allow gtest to pretty print values of type 'dna_base':
inline void PrintTo(const dna_base& base, ::std::ostream* os) {
	(*os) << char(base);
}

// Enable expansion of test sequences when writing dna_sequences to a std::ostream.
void enable_test_sequence_expansion();

// Disable expansion of test sequences when writing dna_sequences to a
// std::ostream.
void disable_test_sequence_expansion();

// Facility for annotating a sequence of bases in the debug output.
void add_print_seq_annotation(const dna_sequence& seq, const std::string& description);
void add_print_seq_annotation(const std::string& seq_str, const std::string& description);

// Enable printing of annotations sequences when outputting DNA sequences.
void enable_annotated_sequences();


// This matches a dna sequence of a concrete size.
class dna_sequence_matcher {
 public:
  dna_sequence_matcher() = delete;
  dna_sequence_matcher(const dna_sequence_matcher&) = default;
  dna_sequence_matcher(const char* seq_str) : dna_sequence_matcher(std::string(seq_str)) {}
  dna_sequence_matcher(const std::string& seq_str)
      : dna_sequence_matcher(seq_str, seq_str.size()) {}
  dna_sequence_matcher(const std::string& regex, size_t size) : m_regex(regex), m_size(size) {
    try {
      m_simple_seq.emplace(regex);
    } catch(const io_exception&) {
    }
  }

  testing::Matcher<dna_sequence> matcher() const;

  const boost::optional<dna_sequence>& get_simple() const { return m_simple_seq; }

  size_t size() const { return m_size; }
  friend std::ostream& operator<<(std::ostream& os, const dna_sequence_matcher& m) {
    if (m.m_simple_seq) {
      return os << m.m_simple_seq;
    } else {
      return os << "[sequence of size " << m.m_size << " matching " << m.m_regex << "]";
    }
  }

 private:
  boost::optional<dna_sequence> m_simple_seq;
  const boost::regex m_regex;
  const size_t m_size;
};

namespace dna_testutil {

// Short aliases to generate test sequences.
inline dna_sequence tseq(const std::string& str) {
	return dna_test_sequence(str);
}

inline dna_sequence tseq_rc(const std::string& str) {
	return dna_test_sequence(str).rev_comp();
}

extern dna_sequence dna_A;
extern dna_sequence dna_C;
extern dna_sequence dna_G;
extern dna_sequence dna_T;

dna_sequence drop_front(unsigned num_bases, const dna_sequence& seq);
dna_sequence drop_back(unsigned num_bases, const dna_sequence& seq);

// This matcher matches a dna sequence with the given bases.
dna_sequence_matcher SequenceIs(const std::string& seq_str);

// This matcher matches a long dna sequence of length len that matches
// the regular expression in seq_str.
dna_sequence_matcher LongSequenceMatches(const std::string& seq_regex, size_t len);

}  // namespace
