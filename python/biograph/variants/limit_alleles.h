#pragma once

#include <pybind11/pybind11.h>

#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

void bind_limit_alleles(pybind11::module&);
