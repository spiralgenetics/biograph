#include "modules/io/uuid.h"
#include <boost/regex.hpp>
#include <gtest/gtest.h>

static const boost::regex uuid_regex(
    "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$");

TEST(uuid_test, uuid_test) {
  std::map<char, int> got_char;
  std::set<std::string> distinct_uuids;
  for (int i=0; i < 1000; ++i) {
    auto uuid = make_uuid();
    EXPECT_TRUE(boost::regex_match(uuid, uuid_regex));
    for (char c : uuid) {
      if (c != '-') {
        got_char[c]++;
      }
    }
    distinct_uuids.insert(uuid);
  }

  EXPECT_EQ(1000, distinct_uuids.size());

  // Make sure each hex digit is actually used.
  EXPECT_EQ(got_char.size(), 16);
}
