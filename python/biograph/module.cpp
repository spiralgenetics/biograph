#include <pybind11/pybind11.h>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/ref.hpp>
#include <iostream>
#include <stdexcept>
#include <string>

//#include "python/biograph/internal/module.h"
#include "python/biograph/biograph.h"
#include "python/biograph/biograph_metadata.h"
#include "python/biograph/log.h"
#include "python/biograph/readmap.h"
#include "python/biograph/seqset.h"
#include "modules/io/track_mem.h"

using namespace pybind11;

void bind_dna_sequence(module& m);
void bind_reference(module& m);
void bind_version(module& m);
void bind_variants_module(module& m);
void bind_internal_module(module& m);

#define XCONCAT3(X, Y, Z) X ## Y ## Z
#define CONCAT3(X, Y, Z) XCONCAT3(X, Y, Z)
#define CAPI_NAME(CAPI) CONCAT3(CAPI, PY_MAJOR_VERSION, PY_MINOR_VERSION)
PYBIND11_MODULE(CONCAT3(_capi_, PY_MAJOR_VERSION, PY_MINOR_VERSION), m) {
  // pybind11::handle<> error(PyObject_CallFunctionObjArgs(std::exception,
  // NULL));
  // docstring_options local_docstring_options(true, false, false);
  // This still doesn't work.
  // http://yangacer.blogspot.com/2014/06/exception-translation-from-c-to-python.html
  // ugh
  // https://mail.python.org/pipermail/cplusplus-sig/2012-June/016653.html
  // wrap::exception<std::runtime_error>("RuntimeError", init<std::string>())
  //.def("__str__", &std::runtime_error::what)
  //;

  m.attr("__doc__") = "Internal C API calls for BioGraph SDK";

  bind_readmap(m);
  bind_seqset(m);
  bind_reference(m);
  bind_dna_sequence(m);
  bind_version(m);
  bind_biograph(m);
  bind_biograph_metadata(m);
  bind_logging(m);

  bind_variants_module(m);
  bind_internal_module(m);

  // Limit memory used through python API to 1 petabyte by default.
  set_maximum_mem_bytes(1000ULL /* K */ * 1000ULL /* M */ * 1000ULL /* G */ * 1000ULL /* T */ *
                        1000ULL /* P */);
}
