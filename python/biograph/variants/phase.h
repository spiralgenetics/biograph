#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

void bind_phases(pybind11::module&);
