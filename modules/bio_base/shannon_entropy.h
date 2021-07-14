#pragma once

#include <array>
#include <boost/optional.hpp>
#include <deque>

#include "modules/bio_base/dna_base.h"
#include "modules/bio_base/dna_sequence.h"

// shannon_entropy allows us to determine the number of bases of a
// sequence that we need to reach a minimum entropy threshold.
//
// The entropy threshold provided in a scale that approximates a
// number of bases.  So, for a random sequence, the length and entropy
// will be very close.
//
// Specifying a threshold of, say,
// std::numeric_limits<unsigned>::max() would allow one to calculate
// the entropy of a specific sequence, as long as the sequence is
// smaller than k_max_size.
//
// Implementation:
//
// We calculate the shannon entropy, or the approximate number of bits
// of entropy per kmer of size k_kmer_size, using:
//
// entropy_bits_per_kmer = sum(- p log(p))
//
// over all symbols (where a symbol is a kmer of size k_kmer_size),
// where p is the probability (occurances of symbol divided by
// total symbols seen) of each symbol.  unused symbols are ignored.
//
// We also want to calculate total_entropy = entropy_bits_per_kmer *
// tot_symbols.
//
// To expand this formula, we use p = symbol_count / tot_symbols and
// get the following (all logs are base 2):
//
// bits_of_entropy = sum(- (symbol_count/tot_symbols)
// log(symbol_count/tot_symbols))
// = sum(- (symbol_count/tot_symbols) (log(symbol_count) - log(tot_symbols))
// = sum((symbol_count/tot_symbols) (log(tot_symbols) - log(symbol_count))

// Multiplying by tot_symbols, we get:
// total_entropy = sum(symbol_count) (log(tot_symbols) - log(symbol_count)))
// = sum(symbol_count) log(tot_symbols) - sum(symbol_count log(symbol_count))
//
// Since sum(symbol_count) is tot_symbols, we get:
// = tot_symbols log(tot_symbols) - sum(symbol_count log(symbol_count))
//
// We keep sum(symbol_count log(symbol_count)) * k_precision_factor
// in m_sum_countlogcount and update it whenever we add or remove a
// base.

class shannon_entropy {
 public:
  // The entropy threshold given is compared to the number of bases
  // multiplied by the per-base shannon entropy.
  shannon_entropy(unsigned entropy_threshold);
  ~shannon_entropy();

  // Add the given base or bases to the entropy calculator.
  void push_front(dna_base b);
  void push_front(const dna_sequence& b);

  // Returns the number of bases needed to reach the entropy
  // threshold, or boost::none if it's impossible to reach within
  // k_max_size bases.
  boost::optional<unsigned> length_needed() const;

  // Returns the calculated entropy value of the most recent
  // length_needed() bases.
  unsigned calc_entropy() const;

 private:
  // Number of bases to use for a symbol.
  static constexpr unsigned k_kmer_size = 3;
  // Maximum size of sequence.
  static constexpr unsigned k_max_size = 255;
  // Size of symbol table.
  static constexpr unsigned k_num_symbols = (1 << (2 * k_kmer_size));
  // How to represent a "1" using integer math.
  static constexpr uint64_t k_precision_factor =
      (1ULL << 63) / ((k_max_size + 1) * (k_max_size + 1));
  // Marker for length needed that we're unable to fulfill the requirement.
  static constexpr unsigned k_inf_length_needed =
      std::numeric_limits<unsigned>::max();

  static_assert(k_precision_factor <
                    ((~0ULL) / ((k_max_size + 1) * (k_max_size + 1))),
                "too large k_precision_factor would cause overflow");

  // Returns a new symbol id with the given base rotated in.
  static unsigned push_symbol_id(unsigned symbol_id, dna_base b);
  void add_to_count(unsigned n);
  void remove_from_count(unsigned n);

  unsigned m_entropy_threshold = 0;
  unsigned m_length_needed = k_inf_length_needed;
  std::deque<dna_base> m_bases;

  // Symbol ID of the last symbol pushed on the front
  unsigned m_front_symbol_id = 0;
  // Symbol ID of the last symbol popped off the back.
  unsigned m_back_symbol_id = 0;

  // Number of instances of each symbol, indexed by symbol id.
  std::array<uint8_t, k_num_symbols> m_symbol_counts{};

  // Precomputed log base 2 of 8-bit numbers
  static std::array<uint64_t, 256> m_log2;

  // Sum(count log2(count)
  uint64_t m_sum_countlogcount = 0;
  unsigned m_tot_symbol_count = 0;
};
