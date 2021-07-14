#pragma once

// message digest computation from OpenSSL

#include <boost/filesystem.hpp>
#include "base/base.h"

// Compute the digest of a file or string with any method supported by OpenSSL
std::string mdsum(const std::string& in_str, const std::string& method);
std::string mdsum(const boost::filesystem::path& in_file, const std::string& method);

// Helper for running a sha1 digest
std::string sha1sum(const std::string& in_str);
std::string sha1sum(const boost::filesystem::path& in_file);
