
#include <string.h>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <string>
#include <iostream>
#include <sstream>

std::string hexdump(const std::string& input)
{
	std::istringstream is(input);
	std::stringstream ss;
	size_t address = 0;

	ss << std::hex << std::setfill('0');
	while (is.good()) {
		int nread;
		char buf[16];

		for (nread = 0; nread < 16 && is.get(buf[nread]); nread++);
		if (nread == 0) {
			break;
		}

		// Show the address
		ss << std::setw(8) << address;

		// Show the hex codes
		for (int i = 0; i < 16; i++) {
			if (i % 8 == 0) {
				ss << ' ';
			}
			if (i < nread) {
				ss << ' ' << std::setw(2) << ((unsigned)buf[i] & 0xff);
			}
			else {
				ss << "   ";
			}
		}

		// Show printable characters
		ss << "  ";
		for (int i = 0; i < nread; i++) {
			if (buf[i] < 32) {
				ss << '.';
			}
			else {
				ss << buf[i];
			}
		}

		ss << "\n";
		address += 16;
	}
	return ss.str();
}

