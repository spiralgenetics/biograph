#pragma once

#include "modules/mapred/manifest.h"
#include "modules/mapred/path.h"

#include "modules/io/log.h"
#include "modules/io/mem_io.h"
#include "modules/io/make_unique.h"
#include "modules/io/uuid.h"

// holder_concept describe the interface supported by base_chunker
//
// client -> base_chunker::write(k,v)
//		-> holder_concept::write(k,v)
// ...
// n writes happen
// ...
//
// client -> base_chunker::write(k,v)
// 		current holder is too big, let's chunk it!
// 		holder.prep_read() // do final transformation on data
// 		read holder into encoding buffer
// 		write buffer to file_info
// 		holder.clear() // empty the holder
//
// hold_concept cannot be a class that derives from kv_writer.
// The reason is that prep_read is going to modify the underlying data
// so the encoding to the writer cannot be done on the fly (i.e. as soon as ::write(k,v) returns)

class holder_concept
{
public:
	holder_concept(const std::string& serialized_params);
	void update_split(const std::string& key);
	bool oversized(size_t goal_size) const;
	bool legal_split(const std::string& key) const;
	bool split_now(const std::string& key) const;
	void write(const std::string& key, const std::string& value);
	size_t get_num_records() const;
	void prep_read(); //potentially modify the underlying data before read ops start
	void set_file_info(file_info& fi) const;
	size_t get_size() const;
	void clear();
};

template<class holder>
class base_chunker : public kv_sink
{
public:
	base_chunker(
		const std::string& param,
		const path& path,
		const std::string& name_prefix,
		size_t goal_size,
		size_t partition,
		manifest& out,
		const std::string& encoding
	);

	void write(const std::string& key, const std::string& value) override;
	void flush();
	void close() override;
	void allow_split() { m_allow_split = true; }

private:
	void start_chunk();
	void end_chunk();

	path m_root;
	std::string m_name_prefix;
	size_t m_chunk_id;
	size_t m_goal_size;
	size_t m_partition;
	std::unique_ptr<waiter> m_writer;
	holder m_hold1;
	holder m_hold2;
	holder* m_cursort;
	holder* m_lastsort;
	std::unique_ptr<mem_io> m_curread;
	std::unique_ptr<mem_io> m_lastread;
	path m_path;
	manifest& m_out;
	bool m_allow_split;
	std::string m_encoding;
};

template<class holder>
base_chunker<holder>::base_chunker(
		const std::string& param,
		const path& root,
		const std::string& name_prefix,
		size_t goal_size,
		size_t partition,
		manifest& out,
		const std::string& encoding)
	: m_root(root)
	, m_name_prefix(name_prefix)
	, m_chunk_id(0)
	, m_goal_size(goal_size)
	, m_partition(partition)
	, m_hold1(param)
	, m_hold2(param)
	, m_curread(make_unique<mem_io>("",track_alloc("base_chunker:curread")))
	, m_lastread(make_unique<mem_io>("",track_alloc("base_chunker:lastread")))
	, m_out(out)
	, m_allow_split(false)
	, m_encoding(encoding)
{
	m_out.set_encoding(encoding);
	m_cursort = &m_hold1;
	m_lastsort = &m_hold2;
}

template<class holder>
void base_chunker<holder>::write(const std::string& key, const std::string& value)
{
	if ((m_cursort->oversized(m_goal_size) && (m_allow_split || m_cursort->legal_split(key)))
		|| (m_cursort->split_now(key))) {
		m_hold1.update_split(key);
		m_hold2.update_split(key);
		flush();
	}
	m_cursort->write(key, value);
}

template<class holder>
void base_chunker<holder>::flush()
{
	end_chunk();
	start_chunk();
}

template<class holder>
void base_chunker<holder>::close()
{
	flush();
	end_chunk();
}

template<class holder>
void base_chunker<holder>::end_chunk()
{
	if (!m_writer) {
		return;
	}

	std::unique_ptr<waiter> writer(m_writer.release());

	// SPLOG_P(LOG_DEBUG, "base_chunker::writer> Waiting for previous chunk to be written to storage");
	writer->wait();
	// SPLOG_P(LOG_DEBUG, "base_chunker::writer> Wait complete");

	file_info fi;
	fi.file = m_path;
	m_lastsort->set_file_info(fi);
	m_out.add(fi, m_partition);

	m_lastread->clear();
}

template<class holder>
void base_chunker<holder>::start_chunk()
{
	if (m_cursort->get_num_records() == 0) {
		return;
	}

	// SPLOG_P(LOG_DEBUG, "base_chunker::start_chunk> Starting new chunk");
	m_cursort->prep_read(); //perform the last transformation in the data collected so far

	// read the current block of data and encode (compress) it into an in-memory buffer (m_curread)
	auto writer = make_encoder(m_encoding, *m_curread);
	kv_writer tmp_kv_writer(*writer);
	kv_copy(*m_cursort, tmp_kv_writer);
	tmp_kv_writer.close();
	writer->close();
    // SPLOG_P(LOG_DEBUG, "base_chunker::start_chunk> Done encoding chunk. original size: %lu bytes, actual size: %lu bytes",
    // m_cursort->get_size(), m_curread->size()
    // );

	auto name = printstring("%s_%ld", m_name_prefix.c_str(), m_chunk_id++);
	m_path = m_root.append(printstring("%03ld_%s_%s",
		random() % 1000,
		make_uuid().c_str(),
		name.c_str())
	);

	m_writer = std::move(m_path.write_inverted(*m_curread, m_curread->size()));
	std::swap(m_cursort, m_lastsort);
	std::swap(m_curread, m_lastread);
	m_cursort->clear();
}
