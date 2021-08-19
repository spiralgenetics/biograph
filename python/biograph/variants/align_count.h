#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

// Runs "align_count" pipeline against a set of assemblies.

void bind_align_count(pybind11::module& m);
