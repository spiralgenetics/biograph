#include "python/biograph/log.h"
#include <pybind11/pybind11.h>

#include "modules/io/log.h"
#include "python/common.h"

#include <iostream>

using namespace pybind11;

namespace {

void wrapped_set_spiral_logging_target(object python_logging_target) {
  set_spiral_logging_target([python_logging_target](int level, const std::string& msg) {
    gil_scoped_acquire acquire_gil;
    int py_level;
    switch (level) {
      case LOG_EMERG:
        py_level = 50;  // CRITICAL
        break;
      case LOG_ALERT:
        py_level = 50;  // CRITICAL
        break;
      case LOG_CRIT:
        py_level = 50;  // CRITICAL
        break;
      case LOG_ERR:
        py_level = 40;  // ERROR
        break;
      case LOG_WARNING:
        py_level = 30;  // WARNING
        break;
      case LOG_NOTICE:
        py_level = 20;  // INFO
        break;
      case LOG_INFO:
        py_level = 20;  // INFO
        break;
      case LOG_DEBUG:
        py_level = 10;  // DEBUG
        break;
      default:
        py_level = 20;  // INFO
        break;
    }
    try {
      python_logging_target(py_level, msg);
    } catch (const error_already_set&) {
      static unsigned fail_count = 0;

      ++fail_count;
      if (!(fail_count & (fail_count - 1))) {
        std::cerr << "Logging spiral message to python failed:\n";
        PyErr_Print();
      }

      std::cerr << msg << "\n";
    }
  });
}

}  // namespace

void bind_logging(module& m) {
  m.def("set_spiral_logging_target", wrapped_set_spiral_logging_target);
  m.def("log_build_stamp", log_build_stamp);
}
