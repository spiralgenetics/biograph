#pragma once

#include <pybind11/pybind11.h>

#include "modules/bio_base/biograph.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/spiral_file.h"
#include "python/biograph/readmap.h"
#include "python/biograph/reference.h"
#include "python/biograph/seqset.h"

void bind_biograph(pybind11::module& m);
