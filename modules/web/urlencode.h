#ifndef __urlencode_h__
#define __urlencode_h__

#include <string>

// URL encode or decode ALL input's characters
std::string urlencode(const std::string& input);
std::string urldecode(const std::string& input);

// unlike url[de,en]code, urlencode_component encodes everything but the '/' character
std::string urlencode_component(const std::string& input);

#endif
