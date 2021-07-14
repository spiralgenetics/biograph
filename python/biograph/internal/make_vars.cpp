#include "modules/bio_format/make_vars.h"
#include <pybind11/pybind11.h>

using namespace pybind11;
void bind_make_vars(module& m) {
  class_<make_vars>(m, "make_vars")
      .def(init<const std::string&, size_t, size_t, const std::string&, const std::string&>())
      .def("snp", &make_vars::snp)
      .def("random_insert", &make_vars::random_insert)
      .def("repeat_insert", &make_vars::repeat_insert)
      .def("random_delete", &make_vars::random_delete)
      .def("transpose", &make_vars::transpose)
      .def("print_sequence", &make_vars::print_sequence);
}
