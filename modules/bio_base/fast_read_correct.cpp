#include "modules/bio_base/fast_read_correct.h"

namespace {

kmer_t kmer_shift_left(kmer_t orig, unsigned kmer_size, dna_base b) {
  kmer_t result = orig;
  result <<= 2;
  result |= int(b);
  result &= ~(std::numeric_limits<kmer_t>::max() << (kmer_size * 2));
  return result;
}

// Performs read correction down a slice.
// "kmer" is the kmer of the bases up to but not including the first base in
// "input".
void correct_internal(frc_output* result, string_view input, kmer_t kmer, const frc_params& params,
                      unsigned min_good_run_here, unsigned max_corrections,
                      bool require_run_at_end) {
  auto it = input.begin();
  DCHECK(it != input.end());
  frc_kmer kmer_info;
  DCHECK(params.kmer_lookup_f(kmer, &kmer_info));

  if (*it != 'N') {
    kmer_t next_kmer = kmer_shift_left(kmer, params.kmer_size, dna_base(*it));
    while (params.kmer_lookup_f(next_kmer, &kmer_info)) {
      result->corrected.push_back(dna_base(*it));
      result->kmers.push_back(kmer_info);
      ++it;
      if (it == input.end()) {
        return;
      }
      kmer = next_kmer;
      if (*it == 'N') {
        break;
      }
      next_kmer = kmer_shift_left(kmer, params.kmer_size, dna_base(*it));
    }
  }

  if (result->corrected.size() < min_good_run_here) {
    return;
  }

  if (max_corrections == 0) {
    return;
  }

  ++it;

  dna_base_array<frc_output> try_outputs;
  unsigned best_size = 0;
  dna_base best_b;
  size_t input_left = input.end() - it;
  for (dna_base b : dna_bases()) {
    kmer_t try_kmer = kmer_shift_left(kmer, params.kmer_size, b);

    if (!params.kmer_lookup_f(try_kmer, &kmer_info)) {
      continue;
    }

    try_outputs[b].kmers.reserve(input_left);
    try_outputs[b].kmers.push_back(kmer_info);

    if (it != input.end()) {
      correct_internal(&try_outputs[b], string_view(it, input_left), try_kmer, params,
                       params.min_good_run, max_corrections - 1, require_run_at_end);
    }
    if (require_run_at_end && try_outputs[b].corrected.size() < params.min_good_run) {
      continue;
    }

    CHECK_LT(try_outputs[b].corrections, max_corrections);
    if (try_outputs[b].corrected.size() >= best_size) {
      best_size = try_outputs[b].corrected.size() + 1;
      best_b = b;
    }
  }

  if (best_size) {
    result->corrected.push_back(best_b);
    result->corrections++;
    const auto& best_output = try_outputs[best_b];
    result->corrected += best_output.corrected;
    result->corrections += best_output.corrections;
    result->kmers.insert(result->kmers.end(), best_output.kmers.begin(), best_output.kmers.end());
  }

  return;
}

}  // namespace

frc_output fast_read_correct(string_view input, const frc_params& params) {
  frc_output result;
  if (input.size() < params.kmer_size) {
    return frc_output{};
  }

  unsigned max_corrections = params.max_corrections;

  auto it = input.begin();
  kmer_t kmer = 0;
  unsigned initial_kmer_left = params.kmer_size;
  frc_kmer kmer_info;

  while (initial_kmer_left || !params.kmer_lookup_f(kmer, &kmer_info)) {
    if (it == input.end()) {
      // Unsuccessful finding any valid kmers.
      return frc_output{};
    }
    if (*it == 'N') {
      ++it;
      initial_kmer_left = params.kmer_size;
      continue;
    }
    kmer = kmer_shift_left(kmer, params.kmer_size, dna_base(*it));
    ++it;
    if (initial_kmer_left) {
      --initial_kmer_left;
    }
  }

  frc_output right_correct;
  if (it == input.begin() + params.kmer_size) {
    result.corrected = dna_sequence(string_view(input.begin(), params.kmer_size));
    if (it == input.end()) {
      result.kmers.push_back(kmer_info);
      return result;
    } else {
      right_correct.kmers.reserve(input.size() - params.kmer_size);
      right_correct.kmers.push_back(kmer_info);
    }
  } else {
    auto kmer_start = it - params.kmer_size;
    std::string left_to_correct(input.begin(), kmer_start - input.begin());
    dna_sequence left_ok(string_view(kmer_start, params.kmer_size));

    std::reverse(left_to_correct.begin(), left_to_correct.end());
    for (char& c : left_to_correct) {
      if (c != 'N') {
        c = char(dna_base(c).complement());
      }
    }
    frc_output left_correct;
    left_correct.kmers.reserve(input.size() - params.kmer_size);
    correct_internal(&left_correct, left_to_correct, rev_comp(kmer, params.kmer_size), params, 0,
                     max_corrections, false /* don't require a run at the end */);

    if (left_correct.corrected.size() != left_to_correct.size()) {
      // Left correction failed.
      return frc_output{};
    }

    result.corrected = left_correct.corrected.rev_comp();
    result.corrected += left_ok;
    for (auto& ki : left_correct.kmers) {
      ki = ki.as_flipped();
    }
    CHECK(result.kmers.empty());
    result.kmers = std::move(left_correct.kmers);
    std::reverse(result.kmers.begin(), result.kmers.end());
    result.kmers.push_back(kmer_info);
    CHECK_LE(left_correct.corrections, max_corrections);
    max_corrections -= left_correct.corrections;
    result.corrections += left_correct.corrections;
  }

  if (it != input.end()) {
    correct_internal(&right_correct, string_view(it, input.end() - it), kmer, params, 0,
                     max_corrections, true /* require a run at the end */);
  }
  CHECK_LE(right_correct.corrections, max_corrections);
  result.corrected += right_correct.corrected;
  result.corrections += right_correct.corrections;
  if (result.kmers.empty()) {
    result.kmers = std::move(right_correct.kmers);
  } else {
    result.kmers.insert(result.kmers.end(), right_correct.kmers.begin(), right_correct.kmers.end());
  }
  return result;
}
