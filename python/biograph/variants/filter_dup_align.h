#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/filter_dup_align.h"

// Runs "filter_dup_align" pipeline against a set of assemblies.

void bind_filter_dup_align(pybind11::module& m);
