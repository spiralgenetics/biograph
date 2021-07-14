#include <pybind11/pybind11.h>

#include "python/common.h"

using namespace pybind11;

void bind_anchor(module& m);
void bind_assemble(module& m);
void bind_vargraph(module& m);
void bind_make_vars(module& m);

void bind_internal_module(module& m) {
  bind_anchor(m);
  bind_assemble(m);
  bind_vargraph(m);
  // bind_make_vars(m);
}
