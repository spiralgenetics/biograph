#pragma once

#include <string>

struct path_impl;
// If s3_io isn't linked, this will be null.
extern path_impl* (*new_path_s3_impl)(const std::string&);
