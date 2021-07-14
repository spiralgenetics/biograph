#pragma once

#include <pybind11/pybind11.h>
#include "modules/variants/assemble.h"

void bind_variants_module(pybind11::module& m);
