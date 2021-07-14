#include <algorithm>
#include <functional>
#include <string>
#include <tr1/array>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include "gtest/gtest.h"

#include "modules/mapred/ex_im_porter_data.h"
#include "modules/bio_format/fastq.h"
#include "modules/bio_format/sam_type.h"
#include "modules/io/file_io.h"
#include "modules/io/keyvalue.h"
#include "modules/io/utils.h"
#include "modules/bio_format/importer.h"
#include "modules/mapred/kv_hold.h"
#include "modules/mapred/pipe_mapper.h"
#include "modules/mapred/pipe_params.h"
#include "modules/mapred/unix_pipeline.h"
#include "modules/test/local_context.h"
#include "modules/test/test_utils.h"
#include "modules/test/fastq_test_utils.h"

std::string make_random(unsigned int size)
{
	std::string r;
	for(unsigned int i = 0; i < size; i++)
		r.push_back('a' + random() % 26);
	return r;
}

class delim_importer : public importer
{
public:
	delim_importer(readable & source, char delim = ' ')
		: m_source(source), m_delim(delim) {}

  void import(kv_sink& sink, simple_metadata& meta) override;

private:
	readable & m_source;
	char m_delim;
};


void delim_importer::import(kv_sink& sink, simple_metadata& meta)
{
	EXPECT_FALSE(m_delim == '\n');	// The class does not work with newline delimiters.

	const size_t k_read_size = 256;

	typedef std::vector<std::string> TStringVec;
	TStringVec split_data;
	std::vector<char> buffer(k_read_size);
	size_t read_here_index = 0;

	while (size_t read_count = m_source.read(&buffer[read_here_index], k_read_size))
	{
		buffer.resize(read_here_index + read_count);
		split_data.clear();
		boost::split(split_data, buffer, std::bind2nd(std::equal_to<char>(), '\n'));

		std::copy(split_data.rbegin()->begin(), split_data.rbegin()->end(), buffer.begin());
		read_here_index = split_data.rbegin()->size();
		buffer.resize(read_here_index + k_read_size);
		split_data.pop_back();

		for (TStringVec::const_iterator iter = split_data.begin();
			iter != split_data.end();
			++iter)
		{
			// Assumption: neither the key nor the value strings contain the delimeter or a newline.
			std::string::size_type delim_index = iter->find(m_delim);
			EXPECT_FALSE(delim_index == std::string::npos);
			std::string key(iter->begin(), iter->begin() + delim_index);
			std::string value(iter->begin() + delim_index + 1, iter->end());
			sink.write(key, value);
		}
	}
}


class delim_exporter : public exporter
{
public:
    delim_exporter(writable & byte_sink, char delim = ' ')
		: exporter(byte_sink), m_delim(delim) {}

    void write(const std::string & key, const std::string & value) override
    {
        m_sink.write(key.data(), key.size());
        m_sink.write((&m_delim), sizeof(m_delim));
        m_sink.write(value.data(), value.size());
		m_sink.write("\n", 1);
    }

private:
	char m_delim;
};


TEST(PipeMapper, Basic)
{
	const unsigned int k_kv_count = 10000;

	kv_hold	kv_storage("");
	for(unsigned int i = 1; i <= k_kv_count; i++)
	{
		kv_storage.write(boost::lexical_cast<std::string>(i), make_random(i));
	}

	pipe_mapper_buffer the_pipe_mapper_buffer(kv_storage, NULL);
	delim_importer space_delim_importer(the_pipe_mapper_buffer);

	std::vector<std::string> args;
	args.push_back("{print $1, length($2)}");

	unix_pipeline get_value_length(
		the_pipe_mapper_buffer,
		"/usr/bin/awk",
		args);
	delim_exporter a_delim_exporter(get_value_length);
	pipe_mapper the_pipe_mapper(the_pipe_mapper_buffer, a_delim_exporter, space_delim_importer, get_value_length);
	kv_hold kv_processed("");

	the_pipe_mapper(kv_processed);

	SPLOG("Validating results");
	std::string key;
	std::string value;
	for(unsigned int i = 1; i <= k_kv_count; i++)
	{
		kv_processed.read(key, value);
		ASSERT_EQ(key, value);
	}
	SPLOG("%d values verified.", k_kv_count);
}

#ifdef D_CLOUD
// Align some reads with external bwa mem using the pipe mapper
TEST(PipeMapper, Step)
{
	// Full path to bwa binary
	std::string bwa = "/opt/sentieon/libexec/bwa-orig";

	// Import some unaligned reads
	SPLOG("Making yeast kvp from fastq");
	path(make_path("PipeMapper/step/")).mkdir();
	path yeast_kvp_path(make_path("PipeMapper/step/yeast_10000.kvp"));
	make_fastq_kv("golden/ftest/yeast_10000.fq", yeast_kvp_path.bare_path());

	SPLOG("Adding reads");
	manifest yeast_unaligned_reads;
	yeast_unaligned_reads.add(file_info(yeast_kvp_path.bare_path(), 1917708, 10000), 0);

	SPLOG("Setting up exporter");
	ex_im_porter_data exporter_data;
	exporter_data.ref_name = "saccharomyces_cerevisiae_EF4";

	SPLOG("Setting up pipe params");
	pipe_params the_pipe_params;
	the_pipe_params.command = "/opt/spiral/wrappers/test_wrapper.py";
	the_pipe_params.args.push_back(bwa);
	the_pipe_params.args.push_back("mem");
	the_pipe_params.args.push_back("-M");
	the_pipe_params.args.push_back("-t");
	the_pipe_params.args.push_back("16");
	the_pipe_params.args.push_back("datasets/reference/saccharomyces_cerevisiae_EF4/source.fasta");
	the_pipe_params.args.push_back("-");
	the_pipe_params.working_dir = "";
	// This is bwa mem's STDIN
	the_pipe_params.exporter_type = "fastq";
	// This is bwa mem's STDOUT
	the_pipe_params.importer_type = "sam";
	the_pipe_params.ex_im_porter_data = msgpack_serialize(exporter_data);

	SPLOG("Setting up wrapper command:");
	SPLOG("%s %s", the_pipe_params.command.c_str(), boost::join(the_pipe_params.args, " ").c_str());
	SPLOG("Working directory: %s", the_pipe_params.working_dir.c_str());
	SPLOG("Importer: %s", the_pipe_params.importer_type.c_str());
	SPLOG("Exporter: %s", the_pipe_params.exporter_type.c_str());

	// Temporarily override reference_path
	std::string refpath = CONF_S(reference_path);
	Config::set("reference_path", "datasets/reference");

	SPLOG("Setting up context");
	local_context context(2, 1000000, path(make_path("PipeMapper/step/bulkdata")));

	manifest yeast_aligned_reads = context.map_only(
			"pipe", json_serialize(the_pipe_params),
			yeast_unaligned_reads, true);

	std::string output_sam_path = make_path("PipeMapper/step/yeast_aligned_reads.sam");
	file_writer output_sam(output_sam_path);

	manifest_reader the_manifest_reader(yeast_aligned_reads);
	sam_exporter the_sam_exporter(output_sam, "saccharomyces_cerevisiae_EF4", false, "", "");
	the_sam_exporter.export_from(the_manifest_reader);

	the_sam_exporter.close();
	output_sam.close();

	Config::set("reference_path", refpath);

	std::string result = sha1sum(boost::filesystem::path(output_sam_path));
	ASSERT_EQ("5cda284bf743f2547304c709f84e73a108e77066", result);
}
#endif // D_CLOUD
