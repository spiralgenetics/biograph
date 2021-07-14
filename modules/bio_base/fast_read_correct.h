#pragma once

#include "modules/bio_base/dna_base.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"
#include "modules/io/string_view.h"

#include <functional>

struct frc_kmer {
  // True if kmer was flipped to look up in the kmer set.
  bool flipped;

  // Index in the kmer set of the looked up kmer.
  size_t index;

  bool operator==(const frc_kmer& rhs) const {
    return flipped == rhs.flipped && index == rhs.index;
  }
  frc_kmer as_flipped() const {
    frc_kmer result = *this;
    result.flipped = !result.flipped;
    return result;
  }

  friend std::ostream& operator<<(std::ostream& os, const frc_kmer ki) {
    if (ki.flipped) {
      os << "rev-kmer@";
    } else {
      os << "fwd-kmer@";
    }
    os << ki.index;
    return os;
  };
};

// Performs fast read correction on a sequence.
//
// "input" is the sequence on which to perform read correction.
//
// "max_corrections" specifies the maximum number of corrections that
// should be done.
//
// "kmer_size" specifies the number of bases in each kmer.
//
// "lookup_f" is a function that returns true if a kmer is valid.
//
// Returns the longest correct sequence that can be generated, and a
// list of offsets of the corrections.
struct frc_output {
  dna_sequence corrected;
  unsigned corrections = 0;
  std::vector<frc_kmer> kmers;
};

struct frc_params {
  unsigned max_corrections = 2;
  unsigned min_good_run = 2;
  unsigned kmer_size = 30;
  std::function<bool(kmer_t, frc_kmer*)> kmer_lookup_f;
};

frc_output fast_read_correct(string_view input, const frc_params& params);
