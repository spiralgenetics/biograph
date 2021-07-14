#pragma once

#include <pybind11/pybind11.h>
#include "modules/variants/apply_edges.h"

void bind_apply_edges(pybind11::module& m);
