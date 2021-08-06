#pragma once

#include <pybind11/pybind11.h>
#include "modules/variants/apply_graph.h"

void bind_apply_graph(pybind11::module& m);
