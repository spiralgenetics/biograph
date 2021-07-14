#pragma once

// Set of distances in bases, used for calculating distances between variants.

#include "modules/variants/read_set.h"

#include <iterator>

namespace variants {

class dist_set : private read_id_set {
 public:
  dist_set() = default;
  ~dist_set() = default;

  using read_id_set::begin;
  using read_id_set::end;
  using read_id_set::k_mask_bits;
  using read_id_set::to_vector;
  using read_id_set::operator|;
  using read_id_set::empty;
  using read_id_set::operator=;
  using read_id_set::iterator;
  using read_id_set::size;

  dist_set& operator|=(const dist_set& rhs) {
    *this = *this | rhs;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const dist_set& dists);

  // Returns (target - elem) for the elem that's closest to target.
  // If there are two equidistant elements, returns the lesser of two.
  int32_t closest_distance_to(int32_t target) const;

  // Inserts an individual offset.
  void insert(int32_t offset);

  // Returns a new dist set, with to_add added to each element.
  dist_set add_offset(int32_t to_add, int32_t max_val = std::numeric_limits<int32_t>::max(),
                      int32_t max_ideal_val = std::numeric_limits<int32_t>::max()) const;

  // Merges this dist set with another and returns the result.
  void insert_dists(const dist_set& rhs);

  int32_t max_value() const { return closest_distance_to(std::numeric_limits<int32_t>::max()); }

  int32_t min_value() const { return closest_distance_to(0); }

  bool contains(int32_t dist) const;

 private:
  using mask_t = read_id_mask_t;

  static uint32_t max_val_in(uint32_t chunk_id, mask_t bits);
  static uint32_t min_val_in(uint32_t chunk_id, mask_t bits);
};

}  // namespace variants
