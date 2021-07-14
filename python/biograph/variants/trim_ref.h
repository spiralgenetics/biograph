#pragma once

#include <pybind11/pybind11.h>

#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"
#include "modules/variants/assemble.h"
#include "modules/variants/trim_ref.h"

void bind_trim_ref(pybind11::module& m);
