#include "modules/io/uuid.h"
#include <random>
#include <string>

std::string make_uuid()
{
  // Explicitly specify /dev/urandom due to
  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94087
  std::random_device random_dev("/dev/urandom");

  static const char hex_digits[] = "0123456789abcdef";
  std::uniform_int_distribution<int> hex_digit(0, 15);

  // Generate strings like:
  // 1b4e28ba-2fa1-11d2-883f-0016d3cca427
  std::string result;
  for (auto section_length : {8, 4, 4, 4, 12}) {
    if (!result.empty()) {
      result += "-";
    }
    for (int i=0; i < section_length; ++i) {
      result += hex_digits[hex_digit(random_dev)];
    }
  }
  return result;
}
