#include "modules/web/urlencode.h"
#include "modules/io/log.h"
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>

std::string urlencode(const std::string& input)
{
	// Implementation adapted from boost/network/protocol/http/impl/message.ipp

	char encode_buf[4];
	std::string result;
	encode_buf[0] = '%';
	result.reserve(input.size());

	// character selection for this algorithm is based on the following url:
	// http://www.blooberry.com/indexdot/html/topics/urlencoding.htm

	for (size_t pos = 0; pos < input.size(); ++pos) {
		switch (input[pos]) {
			default:
				if (input[pos] >= 32 && input[pos] < 127) {
					// character does not need to be escaped
					result += input[pos];
					break;
				}
			// else pass through to next case

			case '$':
			case '&':
			case '+':
			case ',':
			// We need to be smarter about when to encode / (eg. only in a query string)
			case '/':
			case ':':
			case ';':
			case '=':
			case '?':
			case '@':
			case '"':
			case '<':
			case '>':
			case '#':
			case '%':
			case '{':
			case '}':
			case '|':
			case '\\':
			case '^':
			case '~':
			case '[':
			case ']':
			case '`':
			case ' ':
				// the character needs to be encoded
				sprintf(encode_buf + 1, "%02X", input[pos]);
				result += encode_buf;
				break;
		}
	};

	return result;
}

std::string urldecode(const std::string& input)
{
	// Implementation taken from boost/network/protocol/http/impl/message.ipp
	char decode_buf[3];
	std::string result;
	result.reserve(input.size());

	for (size_t pos = 0; pos < input.size(); ++pos) {
		switch (input[pos]) {
			case '+':
				// convert to space character
				result += ' ';
				break;
			case '%':
				// decode hexidecimal value
				if (pos + 2 < input.size()) {
					decode_buf[0] = input[++pos];
					decode_buf[1] = input[++pos];
					decode_buf[2] = '\0';
					result += static_cast<char>(strtol(decode_buf, 0, 16));
				} else {
					// recover from error by not decoding character
					result += '%';
				}
				break;
			default:
				// character does not need to be escaped
				result += input[pos];
		}
	};

	return result;
}

std::string urlencode_component(const std::string& input)
{
	std::vector<std::string> components;
	boost::split(components, input, boost::is_any_of("/"), boost::token_compress_on);
	std::transform( components.begin(), components.end(), components.begin(), urlencode );
	return boost::algorithm::join(components, "/");
}
