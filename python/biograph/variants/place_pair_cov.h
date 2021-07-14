#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/place_pair_cov.h"

void bind_place_pair_cov(pybind11::module& m);
