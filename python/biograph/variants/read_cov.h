#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

// Runs "read_cov" pipeline against a set of assemblies.

void bind_read_cov(pybind11::module& m);
