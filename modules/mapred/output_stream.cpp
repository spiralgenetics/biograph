#include "modules/mapred/output_stream.h"

#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/base_chunker.h"
#include "modules/mapred/path.h"
#include "modules/mapred/task.h"
#include "modules/mapred/kv_hold.h"
#include "modules/mapred/kv_summarize.h"
#include "modules/mapred/kv_sort.h"
#include <vector>

template<class holder>
class base_output_stream : public kv_sink
{
public:
	base_output_stream(
		const path& base_path,
		const std::string& name_prefix,
		const output_stream_params& params,
		const std::string& holder_param, manifest& out,
		const std::string& encoding,
		const std::string& begin_on,
		const std::string& end_before,
		bool clean_break
	);

	~base_output_stream();

	// For use as kv_sink
	void write(const std::string& key, const std::string& value) override;
	void close() override;

	void allow_split()
	{
		for (size_t i = 0; i < m_parts.size(); i++) {
			m_parts[i]->allow_split();
		}
	}

private:
	void try_run();
	void cleanup();

private:
	manifest& m_output;
	std::unique_ptr<sorter> m_sorter;
	std::vector<base_chunker<holder>*> m_parts;
	std::string m_begin_on;
	std::string m_end_before;
	bool m_clean_break;
};


template<class holder>
base_output_stream<holder>::base_output_stream(
		const path& base_path,
		const std::string& name_prefix,
		const output_stream_params& params,
		const std::string& holder_param,
		manifest& out,
		const std::string& encoding,
		const std::string& begin_on,
		const std::string& end_before,
		bool clean_break)
	: m_output(out)
	, m_begin_on(begin_on)
	, m_end_before(end_before)
	, m_clean_break(clean_break)
{
	// SPLOG_P(LOG_DEBUG, "Output Stream: Constructing");

	// Clear the manifest
	m_output = manifest(params.sort, params.num_partitions);

	// Set up the sorter
	std::string sort = (params.sort == "") ? "lexical" : params.sort;
	m_sorter = sorter_registry::get(sort, "");
	if (!m_sorter) {
		throw io_exception("Unknown sorter");
	}

	// Set up the 'output' stuff
	// SPLOG_P(LOG_DEBUG, "Output Stream: Opening chunkers");
	for (size_t i = 0; i < params.num_partitions; i++) {
		m_parts.push_back(
			new base_chunker<holder>(
				holder_param,
				base_path,
				name_prefix,
				params.goal_size,
				i,
				m_output,
				encoding
			)
		);
	}
}

template<class holder>
base_output_stream<holder>::~base_output_stream()
{
	for (size_t i = 0; i < m_parts.size(); i++) {
		delete m_parts[i];
	}
}


template<class holder>
void base_output_stream<holder>::write(const std::string& key, const std::string& value)
{
	if (m_begin_on != "" && !(m_clean_break ?
		m_sorter->compare(m_begin_on, key) < 2 :
		m_sorter->compare(m_begin_on, key) < 1)) {
		return;
	}

	if (m_end_before != "" && (m_clean_break ?
		m_sorter->compare(m_end_before, key) < 2 :
		m_sorter->compare(m_end_before, key) < 1)) {
		return;
	}

	size_t part = (m_parts.size() > 1) ? m_sorter->partition(key, m_parts.size()) : 0;
	m_parts[part]->write(key, value);
}

template<class holder>
void base_output_stream<holder>::close()
{
	// DO NOT DO SYNCED CLOSE!
	//SPLOG_P(LOG_DEBUG, "Output Stream: Flushing");
	//for(size_t i = 0; i < m_parts.size(); i++)
	//	m_parts[i]->flush();

	// SPLOG_P(LOG_DEBUG, "Output Stream: Closing");
	for (size_t i = 0; i < m_parts.size(); i++) {
		m_parts[i]->close();
	}

	// SPLOG_P(LOG_DEBUG, "Output Stream: Closed");
}

std::unique_ptr<kv_sink> output_stream_params::build(
	const path& base_path,
	const std::string& name_prefix,
	manifest& out)
{
	if (encoding.empty()) {
		encoding = codec::gzip;
		// SPLOG_P(LOG_DEBUG, "output_stream_params::build> detected empty encoding, forcing encoding to %s", encoding.c_str());
	}
	if (num_partitions == 0) {
		throw io_exception("Must have at least 1 partitions");
	}

	if (sort == "") {
		if (reduce != "") {
			throw io_exception("Can't have a reducer without a sorter");
		}
		if (num_partitions != 1) {
			throw io_exception("num_partitions is invalid for unsorted outputs");
		}
		return make_unique<base_output_stream<kv_hold>>(
			base_path, name_prefix, *this, "", out, encoding, begin_on, end_before, clean_break);
	}

	if (reduce != "") {
		kv_summarize::param p;
		p.sort = sort;
		p.reduce = reduce;
		p.reduce_param = reduce_param;
		return make_unique<base_output_stream<kv_summarize>>(
			base_path, name_prefix, *this, json_serialize(p), out, encoding, begin_on, end_before, clean_break);
	}

	if (presorted) {
		auto ret = make_unique<base_output_stream<kv_hold>>(
			base_path, name_prefix, *this, sort, out, encoding, begin_on, end_before, clean_break);
		if (allow_split) {
			ret->allow_split();
		}
		return std::move(ret);
	}

	if (! split.empty() && begin_on.empty()) {
		SPLOG("Splitters (%s) require that the first key be passed in output_stream_params::begin_on", split.c_str());
		throw io_exception("Splitter has no initial key set in output stream parameters.");
	}

	kv_sort::param sort_param;
	sort_param.sorter = sort;
	sort_param.splitter = split;
	sort_param.first_key = begin_on;

	return make_unique<base_output_stream<kv_sort>>(
		base_path, name_prefix, *this, msgpack_serialize(sort_param), out, encoding, begin_on, end_before, clean_break);
}

