#pragma once

#include <pybind11/pybind11.h>
#include "python/biograph/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/spiral_file.h"

void bind_seqset(pybind11::module& m);
