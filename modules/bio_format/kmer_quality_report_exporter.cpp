#include "modules/bio_format/kmer_quality_report_exporter.h"
#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/io/file_io.h"
#include "modules/io/config.h"
#include "modules/io/json_transfer.h"
#include "modules/io/msgpack_transfer.h"
#include <sstream>
#include <fstream>
#include <inttypes.h>

struct kmer_quality_metadata
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(name, TF_STRICT); 
	}
	std::string name;
};

struct kmer_quality_data
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(metadata, TF_STRICT); 
		FIELD(data, TF_STRICT); 
	}
	kmer_quality_metadata metadata;
	std::vector<point> data;
};

kmer_quality_report_exporter::kmer_quality_report_exporter(
	writable& byte_sink, 
	const std::string& dataset_name) 
	: exporter(byte_sink)
	, m_dataset_name(dataset_name) 
{}

void kmer_quality_report_exporter::write_header()
{
	std::string data_root(CONF_S(install_root) + "/etc/kmer_quality/");
	
	std::vector<std::string> files = {
		"header.html"
	};

	for (auto file : files) {
		file_reader fr(data_root + file);
		io_copy(fr, m_sink);
	}
}

void kmer_quality_report_exporter::write(const std::string& key, const std::string& value)
{
	uint64_t k, v;

	msgpack_deserialize(k, key);
	msgpack_deserialize(v, value);
	point p {k,v};
	m_data.push_back(p);
}

void kmer_quality_report_exporter::write_footer()
{
	std::string assign_str("kmer_data = ");
	m_sink.write(assign_str.data(), assign_str.size());
	
	std::vector<kmer_quality_data> all_data;
	kmer_quality_data kmer_data;
	kmer_data.metadata.name = m_dataset_name.c_str();
	kmer_data.data = m_data; 
	all_data.push_back(kmer_data);
	std::string all_data_str = json_serialize(all_data);
	m_sink.write(all_data_str.data(), all_data_str.size());
	
	std::string data_root(CONF_S(install_root) + "/etc/kmer_quality/");
	
	std::vector<std::string> files = {
		"footer.html"
	};
   
	for (auto file : files) {
		file_reader fr(data_root + file);
		io_copy(fr, m_sink);
	}
}
