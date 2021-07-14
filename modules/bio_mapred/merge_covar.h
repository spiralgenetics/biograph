#pragma once

#include <map>
#include <memory>
#include <cmath>

#include "modules/io/int_seq.h"
#include "modules/io/io.h"

void do_merge_covar(const std::vector<std::string>& m_input_paths, const std::string& m_output_path);

static inline double qual_to_prob(double qual) { return std::pow(10.0, -qual / 10.0); }
static inline double prob_to_qual(double prob) { return -10 * log10(prob); }
static inline double err_to_qual(double errors, uint64_t observations) { return prob_to_qual((errors + 1) / (observations + 2)); }

const std::vector<char> k_event_type_lookup = {'M', 'I', 'D'};

static inline std::string four_sig_fig(double x)
{
	std::array<char, 32> buffer;
	std::snprintf(buffer.data(), buffer.size(), "%.4f", x);
	return std::string(buffer.data());
}

static inline std::string two_sig_fig(double x)
{
	std::array<char, 32> buffer;
	std::snprintf(buffer.data(), buffer.size(), "%.2f", x);
	return std::string(buffer.data());
}

const char k_rec_separator = '\036'; // ASCII Record Separator.

class readable;
class writable;
class abstract_table;
class merge_covar_data
{
public:
	void merge_one_file(const std::string& file_path);
	void output(const std::string& output_path) const;

private:
	std::string m_current_file_path;
	std::vector<std::unique_ptr<abstract_table>> m_tables;
	
	void check_table_version(readable& input_reader) const;
	void output_table_version(writable& output) const;
};


class abstract_table
{
public:
	explicit abstract_table(std::string file_path) : m_file_path(std::move(file_path)) {}
  virtual ~abstract_table() = default;
	
	void parse_header(readable& input_reader);
	
	virtual void merge(readable& new_file) = 0;
	virtual void output(writable& table_target) const = 0;

protected:
	int m_rows;
	int m_columns;
	std::string m_name;   // e.g. "Arguments"
	std::string m_title;  // e.g. "Recalibration argument collection values used in this run"
	std::vector<std::string> m_header; // e.g. {"Argument", "Value"}
	std::string m_file_path;

	void check_field_count(
		const std::string& in_what
		, const std::string& the_line
		, int field_count
		, int expected_count) const;
		
	void check_header_field(
		const std::string& which_field
		, int field_count
		, const std::string& the_line
		, const std::string& actual_field
		, const std::string& expected_field) const;

	std::vector<std::string> read_table_header(readable& file_input) const;
	void check_table_header(
		readable& file_input
		, const std::string& table_name
		, const std::vector<std::string>& existing_header) const;

	void parse_table_header(
		readable& file_input
		, const std::string& table_name
		, bool is_first_time
		, const std::vector<std::string>& existing_header);

	std::string header_to_string(const std::vector<std::string>& table_header) const;
	int find_event_index(char event_type) const;
};


// #:SENTIEON_QCAL_TABLE:Arguments:Recalibration argument collection values used in this run
class arguments : public abstract_table
{
public:
	using entry_t = std::tuple <std::string, std::string>;

	explicit arguments(std::string file_path) : abstract_table(std::move(file_path)) {}
	
	void merge(readable& new_file) override;
	void output(writable& table_target) const override;
private:
	std::vector<entry_t> m_data;

	std::string table_row_string(const entry_t& row) const
	{
		return std::get<0>(row) + "\t" + std::get<1>(row);
	}
};

// #:SENTIEON_QCAL_TABLE:Quantized:Quality quantization map
class quantized : public abstract_table
{
public:
	explicit quantized(std::string file_path) : abstract_table(std::move(file_path)) {}
	
	void merge(readable& new_file) override;
	void output(writable& table_target) const override;
private:
	using count_qscore_t = std::pair<int64_t, int32_t>;
	using data_map_t = std::map<int32_t, count_qscore_t>;
	data_map_t m_data;

	std::string table_row_string(const data_map_t::value_type& row) const
	{
		return std::to_string(row.first) + "\t"
			+ std::to_string(std::get<0>(row.second)) + "\t"
			+ std::to_string(std::get<1>(row.second));
	}
};

// Note that all quality calculations are done in probability space and qualities in the
// tables are zero-normalized.  So p = 10^(-Q/10) where p is the probability and Q is the
// normalized quality score.

// The empirical quality is not stored, but calculated from the errors and observations.
// The EstimatedQReported is stored as a product of the probability calculated from the
// quality and the number of observations normalized by the total number of observations.

// #:SENTIEON_QCAL_TABLE:RecalTable0
class recal0 : public abstract_table
{
public:
	explicit recal0(std::string file_path) : abstract_table(std::move(file_path)) {}
	
	void merge(readable& new_file) override;
	void output(writable& table_target) const override;
private:
	struct map_key
	{
		std::string m_read_group;
		int m_event_type;
		
		map_key(std::string read_group, int event_type)
			: m_read_group(std::move(read_group)), m_event_type(event_type)
			{}

		bool operator<(const map_key& rhs) const
		{
			if (m_read_group < rhs.m_read_group) return true;
			if (m_read_group > rhs.m_read_group) return false;
			return m_event_type < rhs.m_event_type;
		}
	};
	
	using recal0_data_t = std::tuple<double, int64_t, double>;
	using data_map_t = std::map<map_key, recal0_data_t>;
	data_map_t m_data;
	
	std::string table_row_string(const data_map_t::value_type& row) const
	{
		auto observations = std::get<1>(row.second);
		return row.first.m_read_group + "\t"
			+ k_event_type_lookup[row.first.m_event_type] + "\t"
			+ two_sig_fig(err_to_qual(std::get<2>(row.second), std::get<1>(row.second))) + "\t"
			+ four_sig_fig(prob_to_qual(std::get<0>(row.second) / observations)) + "\t"
			+ std::to_string(observations) + "\t"
			+ two_sig_fig(std::get<2>(row.second));
	}
};

// #:SENTIEON_QCAL_TABLE:RecalTable1
class recal1 : public abstract_table
{
public:
	explicit recal1(std::string file_path) : abstract_table(std::move(file_path)) {}
	
	void merge(readable& new_file) override;
	void output(writable& table_target) const override;
private:
	struct map_key
	{
		std::string m_read_group;
		int32_t m_quality_score;
		int32_t m_event_type;
		
		map_key(std::string read_group, int quality_score, int event_type)
			: m_read_group(std::move(read_group)), m_quality_score(quality_score), m_event_type(event_type)
			{}

		bool operator<(const map_key& rhs) const
		{
			if (m_read_group < rhs.m_read_group) return true;
			if (m_read_group > rhs.m_read_group) return false;
			if (m_quality_score < rhs.m_quality_score) return true;
			if (m_quality_score > rhs.m_quality_score) return false;
			return m_event_type < rhs.m_event_type;
		}
	};
	
	using recal1_data_t = std::tuple<int64_t, double>;
	using data_map_t = std::map<map_key, recal1_data_t>;
	data_map_t m_data;

	std::string table_row_string(const data_map_t::value_type& row) const;
};

// #:SENTIEON_QCAL_TABLE:RecalTable2
class recal2 : public abstract_table
{
public:
	explicit recal2(std::string file_path) : abstract_table(std::move(file_path)) {}
	
	void merge(readable& new_file) override;
	void output(writable& table_target) const override;

	struct map_key
	{
		std::string m_read_group;
		std::string m_covar_value;
		std::string m_covar_name;
		std::string m_source_file;
		int32_t m_quality_score;
		int32_t m_event_type;

		map_key(
			std::string read_group
			, int quality_score
			, std::string covar_value
			, std::string covar_name
			, int event_type
			, std::string file_name)
				: m_read_group(std::move(read_group))
				, m_covar_value(std::move(covar_value))
				, m_covar_name(std::move(covar_name))
				, m_source_file(std::move(file_name))
				, m_quality_score(quality_score)
				, m_event_type(event_type)
			{}

		// This table's sorting is FUBAR. The covariate value is a DNA
		// string for the "Context" covariate and a signed integer. To
		// add insult to injury, the DNA strings are sorted as if they
		// are reversed.  WTF Sentieon?
		
		bool reverse_less(const std::string& lhs, const std::string& rhs) const
		{
			auto lhs_iter = lhs.crbegin();
			auto rhs_iter = rhs.crbegin();
			while (lhs_iter != lhs.crend() && rhs_iter != rhs.crend()) {
				if (*lhs_iter < *rhs_iter) return true;
				if (*lhs_iter++ > *rhs_iter++) return false;
			}
			
			return rhs_iter != rhs.crend();
		}
		
		bool operator<(const map_key& rhs) const
		{
			if (m_read_group < rhs.m_read_group) return true;
			if (m_read_group > rhs.m_read_group) return false;
			if (m_covar_name < rhs.m_covar_name) return true;
			if (m_covar_name > rhs.m_covar_name) return false;
			if (m_quality_score < rhs.m_quality_score) return true;
			if (m_quality_score > rhs.m_quality_score) return false;
			if (m_event_type < rhs.m_event_type) return true;
			if (m_event_type > rhs.m_event_type) return false;
			if (m_covar_name == "Context") {
				if (reverse_less(m_covar_value, rhs.m_covar_value)) return true;
				if (m_covar_value != rhs.m_covar_value) return false;
			} else if (m_covar_name == "Cycle") {
				if (std::stoi(m_covar_value) < std::stoi(rhs.m_covar_value)) return true;
				if (std::stoi(m_covar_value) > std::stoi(rhs.m_covar_value)) return false;
			} else {
				std::string error_string = "recal2::operator< running file \"" + m_source_file;
				error_string += "\". Unexpected covariate name found: \"" + m_covar_name + "\".";
				throw io_exception(error_string);
			}
			
			return false;
		}
	};
	
private:
	using recal2_data_t = std::tuple<int64_t, double>;
	using data_map_t = std::map<map_key, recal2_data_t>;
	data_map_t m_data;

	std::string table_row_string(const data_map_t::value_type& row) const;
};
