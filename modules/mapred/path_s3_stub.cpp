#include "modules/mapred/path_s3_stub.h"

path_impl* (*new_path_s3_impl)(const std::string&) = nullptr;
