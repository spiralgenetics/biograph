#include <vector>
#include <string>

#include <boost/algorithm/string/split.hpp>

#include "base/base.h"
#include "modules/io/log.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/bio_mapred/merge_covar.h"

const std::string k_table_id = "SENTIEON_QCAL_TABLE";
const std::string k_table_version = ".V1:5";
const std::string k_header = "#";
const std::string k_separator = ":";
const size_t k_max_line = 4096;


void do_merge_covar(const std::vector<std::string>& input_paths, const std::string& output_path)
{
	CHECK(! input_paths.empty());
	CHECK(! output_path.empty());
		
	merge_covar_data merged_data;
	
	// TODO: Make this async
	for (const auto& input_path : input_paths) {
		merged_data.merge_one_file(input_path);
	}
	
	merged_data.output(output_path);
}

void merge_covar_data::merge_one_file(const std::string& file_path)
{
	m_current_file_path = file_path;
	file_reader input_file(file_path);
	check_table_version(input_file);
	
	if (m_tables.empty()) {
		m_tables.push_back(make_unique<arguments>(file_path));
		m_tables.push_back(make_unique<quantized>(file_path));
		m_tables.push_back(make_unique<recal0>(file_path));
		m_tables.push_back(make_unique<recal1>(file_path));
		m_tables.push_back(make_unique<recal2>(file_path));
	}
	
	for (const auto& table : m_tables) {
		table->parse_header(input_file);
		table->merge(input_file);
	}
}

void merge_covar_data::check_table_version(readable& input_reader) const
{
	std::string expected_version = k_header + k_separator + k_table_id + k_table_version;
	std::string table_version;
	input_reader.readline(table_version, k_max_line);
	
	if (table_version != expected_version) {
		std::string error_string = "merge_covar_data::check_table_version> "
			"Covariate table version mismatch in file \"";
		error_string += m_current_file_path;
		error_string += "\".  Expected \"";
		error_string += expected_version;
		error_string += "\" but got \"";
		error_string += table_version;
		error_string += "\".";
		throw io_exception(error_string);
	}

	SPLOG("Covariate table version is \"%s\" for file \"%s\"", table_version.c_str(), m_current_file_path.c_str());
}

void merge_covar_data::output(const std::string& output_path) const
{
	file_writer output_writer(output_path);
	output_table_version(output_writer);

	for (const auto& table : m_tables) {
		table->output(output_writer);
	}
}

void merge_covar_data::output_table_version(writable& output) const
{
	output.write(k_header + k_separator + k_table_id + k_table_version + "\n");
}

void abstract_table::parse_header(readable& input_reader)
{
	std::string table_dims;
	input_reader.readline(table_dims, k_max_line);
	
	std::vector<std::string> fields;
	boost::split(fields, table_dims, [](char c) { return c == k_separator[0]; });
	
	check_field_count("header", table_dims, fields.size(), 4);
	check_header_field("header", 0, table_dims, fields[0], k_header);
	check_header_field("header", 1, table_dims, fields[1], k_table_id);
	
	m_columns = std::stoi(fields[2]);
	m_rows = std::stoi(fields[3]);

	std::string table_title;
	input_reader.readline(table_title, k_max_line);
	fields.clear();
	boost::split(fields, table_title, [](char c) { return c == k_separator[0]; });
	
	check_field_count("title", table_title, fields.size(), 4);
	
	check_header_field("title", 0, table_title, fields[0], k_header);
	check_header_field("title", 1, table_title, fields[1], k_table_id);
	
	m_name = fields[2];
	m_title = fields[3];
	
	SPLOG("Parsed header for table \"%s\" entitled \"%s\" in file \"%s\"",
		m_name.c_str(), m_title.c_str(), m_file_path.c_str());
}

void abstract_table::check_field_count(
	const std::string& in_what
	, const std::string& the_line
	, int field_count
	, int expected_count) const
{
	if (field_count != expected_count) {
		std::string error_string = "abstract_table::check_field_count> File \"" + m_file_path;
		error_string += "\". Unexpected " + in_what + " field count of ";
		error_string += std::to_string(field_count);
		error_string += " in \"";
		error_string += the_line;
		error_string += "\".  Expected ";
		error_string += std::to_string(expected_count);
		error_string += " fields.";
		throw io_exception(error_string);
	}
}

void abstract_table::check_header_field(
	const std::string& which_field
	, int field_count
	, const std::string& the_line
	, const std::string& actual_field
	, const std::string& expected_field) const
{
	if (actual_field != expected_field) {
		std::string error_string = "abstract_table::check_header_field> File \"" + m_file_path;
		error_string += "\". Unexpected " + which_field + " field ";
		error_string += std::to_string(field_count);
		error_string += " \"";
		error_string += actual_field;
		error_string += "\" in \"" + the_line + "\".  Expected \"";
		error_string += expected_field + "\".";
		throw io_exception(error_string);
	}
}

void abstract_table::check_table_header(
	readable& file_input
	, const std::string& table_name
	, const std::vector<std::string>& existing_header) const
{
	std::vector<std::string> fields = read_table_header(file_input);

	if (existing_header != fields) {
		std::string error_string = "arguments::merge> File \"" + m_file_path;
		error_string += "\". Unexpected " + table_name + " table header. ";
		error_string += "Expected \"" + header_to_string(existing_header) + "\", but found \"";
		error_string += header_to_string(fields) + "\".";
		throw io_exception(error_string);
	}
}


void abstract_table::parse_table_header(
	readable& file_input
	, const std::string& table_name
	, bool is_first_time
	, const std::vector<std::string>& existing_header)
{
	if (is_first_time) {
		m_header = read_table_header(file_input);
	} else {
		check_table_header(file_input, "arguments", m_header);
	}
}

std::vector<std::string> abstract_table::read_table_header(readable& file_input) const
{
	std::vector<std::string> fields;

	std::string table_header;
	file_input.readline(table_header, k_max_line);

	boost::split(fields, table_header, [](char c) { return c == '\t'; });
	
	return fields;
}

std::string abstract_table::header_to_string(const std::vector<std::string>& table_header) const
{
	std::string header_string;
	for (const auto& header_element : table_header) {
		header_string += header_element + '\t';
	}
	header_string.pop_back();
	return header_string;
}

int abstract_table::find_event_index(char event_type) const
{
	auto event_type_iter = std::find(
		k_event_type_lookup.begin()
		, k_event_type_lookup.end()
		, event_type);
	if (event_type_iter == k_event_type_lookup.end()) {
		std::string error_string = "recal0::merge> While merging the recal0 table in file \"" + m_file_path;
		error_string += "\" unexpected event type \"";
		error_string += event_type;
		error_string += "\" was found.";
		throw io_exception(error_string);
	}
	
	return event_type_iter - k_event_type_lookup.begin();
}

void arguments::merge(readable& file_input)
{
	bool is_first_time = m_data.empty();
	
	parse_table_header(file_input, "arguments", is_first_time, m_header);
	
	std::string line;
	int line_num = 0;
	std::vector<std::string> fields;
	file_input.readline(line, k_max_line);
	while(! line.empty()) {
		fields.clear();
		boost::split(fields, line, [](char c) { return c == '\t'; });
		check_field_count("argument data", line, fields.size(), 2);
		if (is_first_time) {
			m_data.emplace_back(fields[0], fields[1]);
		} else {
			if (std::get<0>(m_data[line_num]) != fields[0] || std::get<1>(m_data[line_num]) != fields[1]) {
				std::string error_string = "arguments::merge> File \"" + m_file_path;
				error_string = "\". Argument table mismatch on line " + std::to_string(line_num);
				error_string += ". Expected \"" + table_row_string(m_data[line_num]);
				error_string += "\", but found \"" + fields[0] + "\t" + fields[1] +"\".";
			}
		}
		
		file_input.readline(line, k_max_line);
		line_num++;
	}
}

void arguments::output(writable& table_target) const
{
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ std::to_string(m_columns) + k_separator + std::to_string(m_rows)
		+ ";\n"); // No idea why that semicolon is there.  A Sentieon bug perhaps?
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ m_name + k_separator + m_title + "\n");
	table_target.write(header_to_string(m_header) + "\n");
	for (const auto& row : m_data) {
		table_target.write(table_row_string(row) + "\n");
	}
	table_target.write("\n");
}

void quantized::merge(readable& file_input)
{
	bool is_first_time = m_data.empty();
	
	parse_table_header(file_input, "quantized", is_first_time, m_header);
	
	std::string line;
	int line_num = 0;
	std::vector<std::string> fields;
	file_input.readline(line, k_max_line);
	while(! line.empty()) {
		fields.clear();
		boost::split(fields, line, [](char c) { return c == '\t'; });
		check_field_count("quantized data", line, fields.size(), 3);

		int32_t quality = std::stoi(fields[0]);
		int64_t count = std::stol(fields[1]);
		int32_t qscore = std::stoi(fields[2]);
		
		if (is_first_time) {
			m_data.emplace_hint(m_data.end(), quality, std::make_pair(count, qscore));
		} else {
			auto current_row = m_data.find(quality);
			if (current_row == m_data.end()) {
				std::string error_string = "While merging quantized table in file \"" + m_file_path;
				error_string += "\" quality " + fields[0] + " was expected but not found.";
				throw io_exception(error_string);
			}
			
			std::get<0>(current_row->second) += count;
			if (std::get<1>(current_row->second) == 93) {
				std::get<1>(current_row->second) = qscore; // TODO: Figure out what to really do here.
			} else {
				// TODO: Figure out what to do here.
			}
		}
		
		file_input.readline(line, k_max_line);
		line_num++;
	}
}

void quantized::output(writable& table_target) const
{
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ std::to_string(m_columns) + k_separator + std::to_string(m_data.size())
		+ "\n");
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ m_name + k_separator + m_title + "\n");
	table_target.write(header_to_string(m_header) + "\n");
	for (const auto& row : m_data) {
		table_target.write(table_row_string(row) + "\n");
	}
	table_target.write("\n");
}

void recal0::merge(readable& file_input)
{
	bool is_first_time = m_data.empty();
	
	parse_table_header(file_input, "recal0", is_first_time, m_header);

	std::string line;
	int line_num = 0;
	std::vector<std::string> fields;
	file_input.readline(line, k_max_line);
	while(! line.empty()) {
		fields.clear();
		boost::split(fields, line, [](char c) { return c == '\t'; });
		check_field_count("recal0 data", line, fields.size(), 6);
		
		const std::string& event_type = fields[1];
		if (event_type.size() != 1) {
			std::string error_string = "recal0::merge> While merging the recal0 table in file \"" + m_file_path;
			error_string += "\" event type \"" + event_type + "\" was found with size ";
			error_string += std::to_string(event_type.size()) + ". Expected a size of 1.";
			throw io_exception(error_string);
		}
		
		map_key key(fields[0], find_event_index(event_type[0]));

		int64_t observations = std::stol(fields[4]);
		double errors = std::stod(fields[5]);
		double estimated = qual_to_prob(std::stod(fields[3])) * observations;

		if (is_first_time) {
			m_data.emplace_hint(m_data.end(), key
				, std::make_tuple(estimated, observations, errors));
		} else {
			auto current_row = m_data.find(key);
			if (current_row == m_data.end()) {
				m_data.emplace_hint(m_data.end(), key
					, std::make_tuple(estimated, observations, errors));
			}
			else {
				std::get<0>(current_row->second) += estimated;
				std::get<1>(current_row->second) += observations;
				std::get<2>(current_row->second) += errors;
			}
		}
		
		file_input.readline(line, k_max_line);
		line_num++;
	}
}

void recal0::output(writable& table_target) const
{
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ std::to_string(m_columns) + k_separator + std::to_string(m_data.size())
		+ "\n");
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ m_name + k_separator + m_title + "\n");
	table_target.write(header_to_string(m_header) + "\n");
	for (const auto& row : m_data) {
		table_target.write(table_row_string(row) + "\n");
	}
	table_target.write("\n");
}

std::string recal1::table_row_string(const data_map_t::value_type& row) const
{
	return row.first.m_read_group + "\t"
		+ std::to_string(row.first.m_quality_score) + "\t"
		+ k_event_type_lookup[row.first.m_event_type] + "\t"
		+ two_sig_fig(prob_to_qual(std::get<1>(row.second) / std::get<0>(row.second))) + "\t"
		+ std::to_string(std::get<0>(row.second)) + "\t"
		+ two_sig_fig(std::get<1>(row.second));
}

void recal1::output(writable& table_target) const
{
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ std::to_string(m_columns) + k_separator + std::to_string(m_data.size())
		+ "\n");
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ m_name + k_separator + m_title + "\n");
	table_target.write(header_to_string(m_header) + "\n");
	for (const auto& row : m_data) {
		table_target.write(table_row_string(row) + "\n");
	}
	table_target.write("\n");
}

void recal1::merge(readable& file_input)
{
	bool is_first_time = m_data.empty();
	
	parse_table_header(file_input, "recal1", is_first_time, m_header);

	std::string line;
	int line_num = 0;
	std::vector<std::string> fields;
	file_input.readline(line, k_max_line);
	while(! line.empty()) {
		fields.clear();
		boost::split(fields, line, [](char c) { return c == '\t'; });
		check_field_count("recal1 data", line, fields.size(), 6);
		
		map_key key(fields[0], std::stoi(fields[1]), find_event_index(fields[2][0]));

		int64_t observations = std::stol(fields[4]);
		double errors = std::stod(fields[5]);

		if (is_first_time) {
			m_data.emplace_hint(m_data.end(), key
				, std::make_tuple(observations, errors));
		} else {
			auto current_row = m_data.find(key);
			if (current_row == m_data.end()) {
				m_data.emplace_hint(m_data.end(), key
					, std::make_tuple(observations, errors));
			}
			else {
				std::get<0>(current_row->second) += observations;
				std::get<1>(current_row->second) += errors;
			}
		}
		
		file_input.readline(line, k_max_line);
		line_num++;
	}
}

std::string recal2::table_row_string(const data_map_t::value_type& row) const
{
	return row.first.m_read_group + "\t"
		+ std::to_string(row.first.m_quality_score) + "\t"
		+ row.first.m_covar_value + "\t"
		+ row.first.m_covar_name + "\t"
		+ k_event_type_lookup[row.first.m_event_type] + "\t"
		+ two_sig_fig(prob_to_qual(std::get<1>(row.second) / std::get<0>(row.second))) + "\t"
		+ std::to_string(std::get<0>(row.second)) + "\t"
		+ two_sig_fig(std::get<1>(row.second));
}

void recal2::output(writable& table_target) const
{
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ std::to_string(m_columns) + k_separator + std::to_string(m_data.size())
		+ "\n");
	table_target.write(k_header + k_separator + k_table_id + k_separator
		+ m_name + k_separator + m_title + "\n");
	table_target.write(header_to_string(m_header) + "\n");
	for (const auto& row : m_data) {
		table_target.write(table_row_string(row) + "\n");
	}
	table_target.write("\n");
}

void recal2::merge(readable& file_input)
{
	bool is_first_time = m_data.empty();
	
	parse_table_header(file_input, "recal2", is_first_time, m_header);

	std::string line;
	int line_num = 0;
	std::vector<std::string> fields;
	file_input.readline(line, k_max_line);
	while(! line.empty()) {
		fields.clear();
		boost::split(fields, line, [](char c) { return c == '\t'; });
		check_field_count("recal2 data", line, fields.size(), 8);
		
		map_key key(fields[0], std::stoi(fields[1]), fields[2] ,fields[3]
			, find_event_index(fields[4][0]), m_file_path);

		int64_t observations = std::stol(fields[6]);
		double errors = std::stod(fields[7]);

		if (is_first_time) {
			m_data.emplace_hint(m_data.end(), key
				, std::make_tuple(observations, errors));
		} else {
			auto current_row = m_data.find(key);
			if (current_row == m_data.end()) {
				m_data.emplace_hint(m_data.end(), key
					, std::make_tuple(observations, errors));
			}
			else {
				std::get<0>(current_row->second) += observations;
				std::get<1>(current_row->second) += errors;
			}
		}
		
		file_input.readline(line, k_max_line);
		line_num++;
	}
}
