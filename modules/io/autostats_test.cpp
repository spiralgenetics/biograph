#include "modules/io/autostats.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sstream>

using namespace testing;

struct test_stats : public autostats_base {
  DECLARE_AUTOSTATS(test_stats,              //
                    ((COUNTER, my_counter))  //
                    ((MAX, my_max)));
};

TEST(autostats, simple) {
  test_stats st;
  st.my_counter++;
  st.my_max.add(1);

  EXPECT_EQ(1, st.my_counter);
  EXPECT_EQ(1, st.my_max.value());
}

TEST(autostats, add_together) {
  test_stats st;
  st.my_counter += 2;
  st.my_max.add(3);

  EXPECT_EQ(2, st.my_counter);
  EXPECT_EQ(3, st.my_max.value());

  test_stats st2;
  st2.my_counter += 7;
  st2.my_max.add(13);

  st += st2;
  EXPECT_EQ(9, st.my_counter);
  EXPECT_EQ(13, st.my_max.value());
}

TEST(autostats, print_to_stream) {
  test_stats st;
  st.my_counter += 10;
  st.my_max.add(2);
  std::stringstream result;
  result << st;

  EXPECT_EQ(result.str(), "Stats: my_counter: 10, my_max: 2");
}

TEST(autostats, omit_blank) {
  test_stats st;
  st.my_counter += 10;
  std::stringstream result;
  result << st;

  EXPECT_EQ(result.str(), "Stats: my_counter: 10");
}

TEST(autostats, empty_to_stream) {
  test_stats st;
  std::stringstream result;
  result << st;

  EXPECT_EQ(result.str(), "Stats: (no stats)");
}

TEST(autostats, value_map) {
  test_stats st;
  st.my_counter += 10;
  st.my_max.add(2);
  std::stringstream result;
  result << st;

  EXPECT_THAT(st.value_map(),
              UnorderedElementsAre(Pair("my_counter", 10), Pair("my_max", 2)));
}
