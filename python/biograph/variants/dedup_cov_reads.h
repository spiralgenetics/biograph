#pragma once

#include <pybind11/pybind11.h>

#include "modules/variants/assemble.h"
#include "modules/variants/dedup_cov_reads.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

void bind_dedup_cov_reads(pybind11::module& m);
