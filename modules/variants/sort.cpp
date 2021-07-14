#include "modules/variants/sort.h"
#include "modules/io/parallel.h"

namespace variants {

void sorter::flush() {
  std::sort(m_queued.begin(), m_queued.end(),
            [this](const assembly_ptr& a, const assembly_ptr& b) { return m_compare_f(*a, *b); });
  for (assembly_ptr& a : m_queued) {
    m_output->add(std::move(a));
  }
  m_queued.clear();
}

}  // namespace variants
