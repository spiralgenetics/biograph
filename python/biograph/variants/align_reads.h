#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "modules/variants/align_reads.h"

// Runs "align_reads" pipeline against a set of assemblies.

void bind_align_reads(pybind11::module& m);
