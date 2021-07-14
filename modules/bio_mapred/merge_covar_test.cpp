#include <boost/algorithm/string/split.hpp>

#include <gtest/gtest.h>

#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/bio_mapred/merge_covar.h"
#include "modules/test/test_utils.h"

constexpr size_t k_max_line = 4096;

TEST(merge_covar, merge)
{
	std::vector<std::string> input_paths = {"golden/hiv_covar_table0", "golden/hiv_covar_table1"};
	std::string output_path = "merged_covar";
	
	do_merge_covar(input_paths, output_path);
}

enum class table_id_t : uint8_t 
{
	undefined = 0, arguments, quantized, recal0, recal1, recal2, overflow
};

table_id_t& operator++(table_id_t &table_id) 
{
	using int_type_t = typename std::underlying_type<table_id_t>::type;
	table_id = static_cast<table_id_t>( static_cast<int_type_t>(table_id) + 1 );

	return table_id;
}

TEST(merge_covar, merge_one)
{
	std::vector<std::string> input_paths = {"golden/hiv_covar_table0"};
	std::string output_path = "one_merged_covar";
	
	do_merge_covar(input_paths, output_path);
	
	ASSERT_EQ(input_paths.size(), 1);
	file_reader input_file(input_paths[0]);
	file_reader merged_file(output_path);
	
	std::string input_line;
	std::string merged_line;
	std::vector<std::string> input_fields;
	std::vector<std::string> merged_fields;
	
	table_id_t table_id = table_id_t::undefined;
	bool reading_header = true;
	
	while (input_file.readline(input_line, k_max_line)) {
		ASSERT_TRUE(merged_file.readline(merged_line, k_max_line));
		
		if (reading_header) {
			if (input_line[0] != '#') {
				reading_header = false;
				ASSERT_FALSE(table_id == table_id_t::overflow);
				++table_id;
				ASSERT_TRUE(input_line == merged_line);
			} else {
				ASSERT_TRUE(input_line == merged_line);
			}
		} else {
			if (input_line.empty()) {
				ASSERT_TRUE(merged_line.empty());
				reading_header = true;
				break;
			}
			
			switch (table_id) {
				case table_id_t::arguments:
					ASSERT_TRUE(input_line == merged_line);
					break;
					
				case table_id_t::quantized:
					ASSERT_TRUE(input_line == merged_line);
					break;
				
				case table_id_t::recal0:
					input_fields.clear();
					merged_fields.clear();
					boost::split(input_fields, input_line, [](char c) { return c == '\t'; });
					boost::split(merged_fields, merged_line, [](char c) { return c == '\t'; });
					
					ASSERT_EQ(input_fields.size(), 6);
					ASSERT_EQ(merged_fields.size(), 6);
					ASSERT_EQ(input_fields[0], merged_fields[0]);
					ASSERT_EQ(input_fields[1], merged_fields[1]);
					ASSERT_EQ(input_fields[3], merged_fields[3]);
					ASSERT_EQ(input_fields[4], merged_fields[4]);
					ASSERT_EQ(input_fields[5], merged_fields[5]);
					break;
				
				case table_id_t::recal1:
					input_fields.clear();
					merged_fields.clear();
					boost::split(input_fields, input_line, [](char c) { return c == '\t'; });
					boost::split(merged_fields, merged_line, [](char c) { return c == '\t'; });
					
					ASSERT_EQ(input_fields.size(), 6);
					ASSERT_EQ(merged_fields.size(), 6);
					ASSERT_EQ(input_fields[0], merged_fields[0]);
					ASSERT_EQ(input_fields[1], merged_fields[1]);
					ASSERT_EQ(input_fields[2], merged_fields[2]);
					ASSERT_EQ(input_fields[4], merged_fields[4]);
					ASSERT_EQ(input_fields[5], merged_fields[5]);
					break;
				
				case table_id_t::recal2:
					input_fields.clear();
					merged_fields.clear();
					boost::split(input_fields, input_line, [](char c) { return c == '\t'; });
					boost::split(merged_fields, merged_line, [](char c) { return c == '\t'; });
					
					ASSERT_EQ(input_fields.size(), 8);
					ASSERT_EQ(merged_fields.size(), 8);
					ASSERT_EQ(input_fields[0], merged_fields[0]);
					ASSERT_EQ(input_fields[1], merged_fields[1]);
					ASSERT_EQ(input_fields[2], merged_fields[2]);
					ASSERT_EQ(input_fields[3], merged_fields[3]);
					ASSERT_EQ(input_fields[4], merged_fields[4]);
					ASSERT_EQ(input_fields[6], merged_fields[6]);
					ASSERT_EQ(input_fields[7], merged_fields[7]);
					break;
				
				case table_id_t::overflow:
					break;
				
				default:
					ASSERT_TRUE(false);
			}
		}
	}
}

TEST(merge_covar, map_key)
{
	recal2::map_key key1("Some_Read_Group", 35, "AA", "Context", 0, "Test File");
	recal2::map_key key2("Some_Read_Group", 35, "AA", "Context", 0, "Test File");
	ASSERT_FALSE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 35, "-123", "Cycle", 0, "Test File");
	key2 = recal2::map_key("Some_Read_Group", 35, "-123", "Cycle", 0, "Test File");
	ASSERT_FALSE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 36, "TTA", "Context", 0, "Test File");
	key2 = recal2::map_key("Tome_Read_Group", 35, "TTT", "Context", 0, "Test File");
	ASSERT_TRUE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 3, "TTA", "Context", 0, "Test File");
	key2 = recal2::map_key("Some_Read_Group", 23, "TTT", "Context", 0, "Test File");
	ASSERT_TRUE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 23, "GCA", "Context", 0, "Test File");
	key2 = recal2::map_key("Some_Read_Group", 23, "AAT", "Context", 0, "Test File");
	ASSERT_TRUE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 35, "-123", "Cycle", 0, "Test File");
	key2 = recal2::map_key("Some_Read_Group", 35, "123", "Cycle", 0, "Test File");
	ASSERT_TRUE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 35, "CA", "Context", 0, "Test File");
	key2 = recal2::map_key("Some_Read_Group", 35, "AC", "Context", 0, "Test File");
	ASSERT_TRUE(key1 < key2);
	ASSERT_FALSE(key2 < key1);

	key1 = recal2::map_key("Some_Read_Group", 35, "AC", "Context", 0, "Test File");
	key2 = recal2::map_key("Some_Read_Group", 35, "AAC", "Context", 0, "Test File");
	ASSERT_TRUE(key1 < key2);
	ASSERT_FALSE(key2 < key1);
}
