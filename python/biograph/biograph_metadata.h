#pragma once

#include <pybind11/pybind11.h>

#include "modules/bio_base/biograph_dir.h"

void bind_biograph_metadata(pybind11::module& m);
