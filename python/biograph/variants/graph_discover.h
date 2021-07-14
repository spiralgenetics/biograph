#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

// Runs "graph_discover" pipeline against a set of assemblies.

void bind_graph_discover(pybind11::module& m);
