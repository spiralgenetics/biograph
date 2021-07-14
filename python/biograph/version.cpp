#include "modules/io/version.h"
#include "tools/build_stamp.h"
#include <pybind11/pybind11.h>

using namespace pybind11;
std::string get_the_sdk_version() {
  return biograph_sdk_current_version.make_string();
}

std::string sdk_get_build_scm_revision() {
  return get_build_scm_revision();
}

void bind_version(module& m) {
  m.def("version", get_the_sdk_version, "BioGraph SDK version.");
  m.def("build_revision", sdk_get_build_scm_revision, "BioGraph build revision.");
}
