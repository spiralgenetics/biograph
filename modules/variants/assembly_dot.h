#pragma once

#include "modules/variants/scaffold.h"

#include <map>
#include <sstream>

namespace variants {

// Generate a dot graph based on assemblies, and show where they link to
// reference.
class assembly_dot {
 public:
  assembly_dot(const scaffold& s);

  void add_assembly(const assembly& a);

  std::string str() {
    if (!m_finalized) {
      finalize();
    }
    return m_result.str();
  }

 private:
  void finalize();
  scaffold m_scaffold;
  int m_id = 0;
  std::multimap<aoffset_t, std::string> m_ref_right_edges;
  std::multimap<aoffset_t, std::string> m_ref_left_edges;
  bool m_finalized = false;

  std::stringstream m_result;
};

}  // namespace variants
