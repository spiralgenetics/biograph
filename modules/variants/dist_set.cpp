#include "modules/variants/dist_set.h"

namespace variants {

int32_t dist_set::closest_distance_to(int32_t signed_target) const {
  uint32_t target = 0;
  if (signed_target > 0) {
    target = signed_target;
  }

  uint32_t chunk_id = target / k_mask_bits;
  unsigned offset = target & (k_mask_bits - 1);

  auto it = std::lower_bound(m_impl.begin(), m_impl.end(), chunk_id);

  if (it == m_impl.end()) {
    CHECK(it != m_impl.begin()) << "Empty dist set?";

    --it;

    return int32_t(max_val_in(it->chunk_id, it->read_id_bits));
  }

  // Hi is the best value that's greater than the target, and
  // Lo is the best value that's less than or equal to the target.
  uint32_t hi_best = std::numeric_limits<uint32_t>::max();
  uint32_t lo_best = std::numeric_limits<uint32_t>::max();

  if (it->chunk_id == chunk_id) {
    // Both greater than and less than could be in this chunk.

    mask_t hi = it->read_id_bits & ((~mask_t(0)) << offset);
    mask_t lo = it->read_id_bits - hi;

    if (lo) {
      lo_best = max_val_in(it->chunk_id, lo);
    } else if (it != m_impl.begin()) {
      auto lo_it = it;
      --lo_it;
      lo_best = max_val_in(lo_it->chunk_id, lo_it->read_id_bits);
    }

    if (hi) {
      hi_best = min_val_in(it->chunk_id, hi);
    } else {
      auto hi_it = it;
      ++hi_it;
      if (hi_it != m_impl.end()) {
        hi_best = min_val_in(hi_it->chunk_id, hi_it->read_id_bits);
      }
    }
  } else {
    if (it != m_impl.begin()) {
      auto lo_it = it;
      --lo_it;
      lo_best = max_val_in(lo_it->chunk_id, lo_it->read_id_bits);
    }
    if (it != m_impl.end()) {
      auto hi_it = it;
      hi_best = min_val_in(hi_it->chunk_id, hi_it->read_id_bits);
    }
  }
  uint32_t best = std::numeric_limits<uint32_t>::max();
  if (lo_best == std::numeric_limits<uint32_t>::max()) {
    CHECK_NE(hi_best, std::numeric_limits<uint32_t>::max());
    best = hi_best;
  } else if (hi_best == std::numeric_limits<uint32_t>::max()) {
    best = lo_best;
  } else {
    CHECK_GE(hi_best, target);
    uint32_t hi_dist = hi_best - target;

    CHECK_LE(lo_best, target);
    uint32_t lo_dist = target - lo_best;

    if (hi_dist < lo_dist) {
      best = hi_best;
    } else {
      best = lo_best;
    }
  }

  CHECK_NE(best, std::numeric_limits<uint32_t>::max());
  return int32_t(best);
}

dist_set dist_set::add_offset(int32_t signed_to_add, int32_t signed_max_val,
                              int32_t signed_max_ideal_val) const {
  CHECK_GE(signed_to_add, 0);
  CHECK_GE(signed_max_val, 0);

  uint32_t to_add = signed_to_add;
  uint32_t max_val = signed_max_val;
  uint32_t max_ideal_val = signed_max_ideal_val;
  uint32_t max_chunk_id = max_val / k_mask_bits;
  uint32_t max_ideal_chunk_id = (max_ideal_val / k_mask_bits) + 1;

  uint32_t to_add_chunk_ids = to_add / k_mask_bits;
  unsigned to_add_offset = to_add & (k_mask_bits - 1);

  mask_t lo_mask = (~mask_t(0)) >> to_add_offset;
  mask_t hi_mask = ~lo_mask;

  dist_set result;
  result.m_impl.reserve(m_impl.size() * (to_add_offset ? 2 : 1));

  for (auto it = m_impl.begin(); it != m_impl.end(); ++it) {
    mask_t lo_bits = it->read_id_bits & lo_mask;
    mask_t hi_bits = it->read_id_bits & hi_mask;

    // lo_bits get shifted into the same chunk id; hi bits get shifted into the next chunk id.
    uint32_t lo_chunk_id = it->chunk_id + to_add_chunk_ids;
    uint32_t hi_chunk_id = it->chunk_id + to_add_chunk_ids + 1;

    mask_t lo_bits_shifted = lo_bits << to_add_offset;
    mask_t hi_bits_shifted = hi_bits >> (k_mask_bits - to_add_offset);
    if (lo_bits) {
      if (lo_chunk_id > max_chunk_id) {
        break;
      }
      CHECK(lo_bits_shifted);
      if (!result.m_impl.empty() && result.m_impl.back().chunk_id == lo_chunk_id) {
        result.m_impl.back().read_id_bits |= lo_bits_shifted;
      } else {
        if (!result.m_impl.empty()) {
          CHECK_LT(result.m_impl.back().chunk_id, lo_chunk_id);
        }
        elem new_elem;
        new_elem.chunk_id = lo_chunk_id;
        new_elem.read_id_bits = lo_bits_shifted;
        result.m_impl.emplace_back(new_elem);
      }
      if (lo_chunk_id > max_ideal_chunk_id) {
        break;
      }
    }

    if (hi_bits) {
      if (hi_chunk_id > max_chunk_id) {
        break;
      }
      CHECK(hi_bits_shifted);
      if (!result.m_impl.empty()) {
        CHECK_LT(result.m_impl.back().chunk_id, hi_chunk_id);
      }

      elem new_elem;
      new_elem.chunk_id = hi_chunk_id;
      new_elem.read_id_bits = hi_bits_shifted;
      CHECK(new_elem.read_id_bits);
      result.m_impl.emplace_back(new_elem);
      if (hi_chunk_id > max_ideal_chunk_id) {
        break;
      }
    }
  }

  return result;
}

void dist_set::insert_dists(const dist_set& rhs) {
  // Let parent class do the work
  (*this) |= rhs;
}

void dist_set::insert(int32_t signed_rhs) {
  CHECK_GE(signed_rhs, 0);

  uint32_t rhs = signed_rhs;
  return ((read_id_set*)this)->insert(rhs);
}

uint32_t dist_set::max_val_in(uint32_t chunk_id, mask_t bits) {
  CHECK_NE(bits, 0);

#define CLZ_VARIANT __builtin_clz
  static_assert(CLZ_VARIANT(~mask_t(0)) == 0,
                "Should use __builtin_clz variant that matches the size of mask_t");
  return chunk_id * k_mask_bits + (k_mask_bits - CLZ_VARIANT(bits) - 1);
#undef CLZ_VARIANT
}

uint32_t dist_set::min_val_in(uint32_t chunk_id, mask_t bits) {
  CHECK_NE(bits, 0);

  return chunk_id * k_mask_bits + __builtin_ctzl(bits);
}

std::ostream& operator<<(std::ostream& os, const dist_set& dists) {
  bool first = true;
  for (int dist : dists) {
    if (first) {
      first = false;
    } else {
      os << ",";
    }
    os << dist;
  }
  if (first) {
    os << "(empty)";
  }
  return os;
}

bool dist_set::contains(int32_t dist) const {
  if (dist < 0) {
    return false;
  }
  return read_id_set::contains(uint32_t(dist));
}

}  // namespace variants
