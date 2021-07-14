#pragma once

#include <vector>
#include <string>

std::string base64_encode(uint8_t const* buf, size_t len);
std::vector<uint8_t> base64_decode(std::string const& str);
