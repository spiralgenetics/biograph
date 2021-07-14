#include "modules/io/autostats.h"

void autostats_base::write_to_stream(std::ostream& os) const {
  bool first = true;
  os << "Stats: ";
  for (const auto& kv : value_map()) {
    if (kv.second == 0) {
      continue;
    }
    if (first) {
      first = false;
    } else {
      os << ", ";
    }
    os << kv.first << ": " << kv.second;
  }
  if (first) {
    os << "(no stats)";
  }
}
