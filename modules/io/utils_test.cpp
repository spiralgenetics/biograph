#include "modules/io/utils.h" 
#include <gtest/gtest.h>

TEST(utils, map_lines)
{
	std::string text = "line1\nline2\nline3\n";
	std::vector<std::string> lines;
	map_lines(text, [&] (const std::string& line) {
		lines.push_back(line);
	});
	ASSERT_EQ(3, lines.size());
	if (lines.size() == 3) {
		ASSERT_EQ("line1", lines[0]);
		ASSERT_EQ("line2", lines[1]);
		ASSERT_EQ("line3", lines[2]);
	}
}

TEST(utils, time_to_RFC3339)
{
	std::tm epoch = { 0 };
	std::string str = time_to_RFC3339(mktime(&epoch));
	ASSERT_EQ("1899-12-31T00:00:00Z", str);
}
