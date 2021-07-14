#pragma once

#include "modules/io/file_io.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_format/struct_var.h"

#include <boost/noncopyable.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stdio.h>


class ipileup;
bool struct_var_adapter(
    const reference& ref,
    const std::function<void(const struct_var& var)>& sink_f,
    const dna_sequence& var_seq, const dna_slice& lower, const dna_slice& upper,
    ipileup* pile, const struct_var& per_assembly, size_t min_depth = 0,
    bool dup_structural = false);

void log_struct_var_adapter_stats();

class make_vars : public boost::noncopyable
{
public:
	make_vars(const std::string& ref_name, size_t leeway, size_t read_length, const std::string& out, const std::string& scaf = "");

	void snp(const std::string& name);
	void random_insert(const std::string& name, size_t size);
	void repeat_insert(const std::string& name, size_t size);
	void random_delete(const std::string& name, size_t size);
	void transpose(const std::string& name, size_t size);
	void print_sequence(FILE* out);

private:
	struct var_info
	{
		std::string name;
		size_t ref_start; // global reference coordinates
		size_t ref_end;   // global reference coordinates
		dna_sequence ref_seq;
		dna_sequence var_seq;
	};

	typedef std::map<size_t, var_info> vars_t;

private:
	size_t random_loc(size_t space = 0);
	void call(const var_info& vi);

private:
	vars_t m_vars;
	reference m_ref;
	boost::mt19937 m_gen;
	size_t m_leeway;
	size_t m_read_length;
	size_t m_next_id;
	file_writer m_writer;
	struct_var_exporter m_exporter;
	std::string m_scaffold;
	size_t m_scaffold_start;
	size_t m_scaffold_len;
};
