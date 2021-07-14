#include "modules/bio_base/shannon_entropy.h"

#include <math.h>
#include <mutex>

constexpr unsigned shannon_entropy::k_kmer_size;
constexpr unsigned shannon_entropy::k_max_size;
constexpr unsigned shannon_entropy::k_num_symbols;
constexpr uint64_t shannon_entropy::k_precision_factor;
constexpr unsigned shannon_entropy::k_inf_length_needed;
std::array<uint64_t, 256> shannon_entropy::m_log2;

shannon_entropy::shannon_entropy(unsigned entropy_threshold)
    : m_entropy_threshold(entropy_threshold) {
  CHECK_GT(entropy_threshold, 1);
  std::once_flag log_table_initialized;
  std::call_once(log_table_initialized, [](){
      for (unsigned i = 0; i < 256; ++i) {
        m_log2[i] = log2(i) * k_precision_factor;
      }
    });
}

shannon_entropy::~shannon_entropy() {
#ifndef NDEBUG
  while (m_bases.size() >= k_kmer_size) {
    m_back_symbol_id =
        push_symbol_id(m_back_symbol_id, *(m_bases.end() - k_kmer_size));
    CHECK_LT(m_back_symbol_id, m_symbol_counts.size());
    CHECK_GT(m_symbol_counts[m_back_symbol_id], 0);
    remove_from_count(m_symbol_counts[m_back_symbol_id]);
    m_symbol_counts[m_back_symbol_id]--;
    add_to_count(m_symbol_counts[m_back_symbol_id]);
    m_bases.pop_back();
  }

  for (const auto& i : m_symbol_counts) {
    CHECK_EQ(i, 0);
  }
#endif
}

void shannon_entropy::push_front(dna_base b) {
  m_bases.push_front(b);
  m_front_symbol_id = push_symbol_id(m_front_symbol_id, b);
  if (m_bases.size() < k_kmer_size) {
    m_back_symbol_id = m_front_symbol_id;
    return;
  }

  CHECK_LT(m_front_symbol_id, m_symbol_counts.size());
  remove_from_count(m_symbol_counts[m_front_symbol_id]);
  m_symbol_counts[m_front_symbol_id]++;
  add_to_count(m_symbol_counts[m_front_symbol_id]);

  if (m_length_needed != k_inf_length_needed) {
    m_length_needed++;
  }

  unsigned entropy = calc_entropy();
  while (entropy >= m_entropy_threshold || m_bases.size() >= k_max_size) {
    CHECK_GT(m_bases.size(), k_kmer_size)
        << "Entropy threshold might be too small: " << m_entropy_threshold;

    if (entropy >= m_entropy_threshold) {
      m_length_needed = m_bases.size();
    } else {
      m_length_needed = k_inf_length_needed;
    }

    m_back_symbol_id =
        push_symbol_id(m_back_symbol_id, *(m_bases.end() - k_kmer_size));
    CHECK_LT(m_back_symbol_id, m_symbol_counts.size());
    CHECK_GT(m_symbol_counts[m_back_symbol_id], 0);
    remove_from_count(m_symbol_counts[m_back_symbol_id]);
    m_symbol_counts[m_back_symbol_id]--;
    add_to_count(m_symbol_counts[m_back_symbol_id]);

    m_bases.pop_back();
    entropy = calc_entropy();
  }
}

unsigned shannon_entropy::push_symbol_id(unsigned symbol_id, dna_base b) {
  symbol_id <<= 2;
  symbol_id += int(b);
  symbol_id &= ~((~0ULL) << (k_kmer_size * 2));
  return symbol_id;
}

void shannon_entropy::remove_from_count(unsigned n) {
  if (!n) {
    return;
  }
  CHECK_GE(m_tot_symbol_count, n);
  m_tot_symbol_count -= n;

  uint64_t countlogcount = m_log2[n] * n;

  CHECK_GE(m_sum_countlogcount, countlogcount);
  m_sum_countlogcount -= countlogcount;
}

void shannon_entropy::add_to_count(unsigned n) {
  if (!n) {
    return;
  }
  m_tot_symbol_count += n;

  uint64_t countlogcount = m_log2[n] * n;

  m_sum_countlogcount += countlogcount;
}

unsigned shannon_entropy::calc_entropy() const {
  if (!m_tot_symbol_count) {
    return 0;
  }

  uint64_t totlogtotcount = m_log2[m_tot_symbol_count] * m_tot_symbol_count;

  CHECK_GE(totlogtotcount, m_sum_countlogcount);

  uint64_t entropy = totlogtotcount - m_sum_countlogcount;
  entropy /= (k_kmer_size * k_precision_factor);
  // 2 bits per base
  entropy >>= 1;

  return entropy;
}

void shannon_entropy::push_front(const dna_sequence& seq) {
  for (auto it = seq.rcbegin(); it != seq.rcend(); ++it) {
    push_front(it->complement());
  }
}

boost::optional<unsigned> shannon_entropy::length_needed() const {
  if (m_length_needed == k_inf_length_needed) {
    return boost::none;
  } else {
    return m_length_needed;
  }
}
