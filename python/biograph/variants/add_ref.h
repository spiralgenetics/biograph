#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

// Runs "add_ref" pipeline against a set of assemblies.

void bind_add_ref(pybind11::module& m);
