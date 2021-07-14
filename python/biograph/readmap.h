#pragma once

#include <pybind11/pybind11.h>
#include "modules/bio_base/seqset_bitmap.h"
#include "modules/bio_base/readmap.h"

void bind_readmap(pybind11::module& m);
