#pragma once

#include "modules/mapred/path.h"
#include "modules/io/log.h"
#include "modules/io/encoding.h"
#include "modules/io/mem_io.h"

#include <future>

// Reads a list of paths sequentially.  DO NOT use both the 'readable' and 'kv_source' methods
// during the same use, since buffering may occur in kv_reader

typedef std::function<void(const path& chunk)> chunk_notify_f;

// The ManifestIterator concept requires the following interface:
// ++ (incrementing)
// * (deferencing) into a path object
// != (inequality)
//
template<class ManifestIterator>
class multi_reader : public readable, public kv_source
{
public:
	multi_reader(
		const ManifestIterator& begin,
		const ManifestIterator& end,
		const std::string& encoding
	);

	// Overrides readable::read
	size_t read(char* buf, size_t len) override;
	// Overrides kv_source::read
	bool read(std::string& key, std::string& value) override;

	void set_notify(const chunk_notify_f& fn)
	{
		m_notify = fn;
	}

private:
	// Advances to the next path; returns false if no more paths
	// are available to read from.
	bool next();

	// Starts prefetches up to read_parallelism at once.
	void start_reads();

private:
	ManifestIterator m_it;
	ManifestIterator m_it_end;
	std::string m_encoding;

	// Maximum number of paths to prefetch at once.
	static const size_t read_parallelism = 16;

	// Maximum size of object to prefetch
	static const size_t maximum_prefetch_size = 128*1024*1024; // 128 MB

	// Paths currently being prefetched.
	struct prefetch_path {
		path m_path;
		std::unique_ptr<readable> m_reader;
	};
	std::deque<std::future<prefetch_path>> m_raw_prefetch;

	// Current raw (encoded) path being read.
	std::unique_ptr<readable> m_raw;
	// Current decoded path being read.
	std::unique_ptr<readable> m_reader;

	// used to transform (char* buf, size_t len) -> (std::string key, std::string value)
	kv_reader m_kv_reader;
	chunk_notify_f m_notify;
};

template<class ManifestIterator>
bool multi_reader<ManifestIterator>::read(std::string& key, std::string& value)
{
	try {
		return m_kv_reader.read(key, value);
	}
	catch (const io_exception& e) {
		throw io_exception(printstring("multi_reader::read> Exception: %s",
			e.message().c_str()
		));
	}
	return false;
}

template<class ManifestIterator>
multi_reader<ManifestIterator>* make_multi_reader(
	const ManifestIterator& begin,
	const ManifestIterator& end,
	const std::string& encoding)
{
	return new multi_reader<ManifestIterator>(begin, end, encoding);
}

template<class ManifestIterator>
multi_reader<ManifestIterator>::multi_reader(
	const ManifestIterator& begin,
	const ManifestIterator& end,
	const std::string& encoding)
	: m_it(begin)
	, m_it_end(end)
	, m_encoding(encoding)
	, m_kv_reader(*this)
{
	start_reads();
	next();
}

template<class ManifestIterator>
void multi_reader<ManifestIterator>::start_reads()
{
	while (m_it != m_it_end && m_raw_prefetch.size() < read_parallelism) {
		path cur_path = *m_it;
		m_it++;

		if (cur_path.size() > maximum_prefetch_size) {
			// Too big to prefetch; wait until we need it to read it.
			std::function<prefetch_path()> reader =
				[cur_path]() {
				return prefetch_path{cur_path, cur_path.read()};
			};
			m_raw_prefetch.push_back(
				std::async(std::launch::deferred, reader));
		} else {
			std::function<prefetch_path()> prefetcher =
				[cur_path]() {
              std::unique_ptr<readable> mem_raw(new mem_io(cur_path.get(), track_alloc("multi_reader:mem_raw")));

				return prefetch_path{cur_path, std::move(mem_raw)};
			};
			m_raw_prefetch.push_back(
				std::async(std::launch::async, prefetcher));
		}
	}
}

template<class ManifestIterator>
bool multi_reader<ManifestIterator>::next()
{

	if (m_raw_prefetch.empty()) {
		// No more parts available!
		return false;
	}

	std::future<prefetch_path> prefetch_future = std::move(m_raw_prefetch.front());
	m_raw_prefetch.pop_front();

	prefetch_path prefetch = prefetch_future.get();

	m_raw = std::move(prefetch.m_reader);
	m_reader = make_decoder(m_encoding, *m_raw);

	// SPLOG_P(LOG_DEBUG, "multi_reader::next> Opening next chunk: %s", prefetch.path.url().c_str());

	if (m_notify) {
		m_notify(prefetch.m_path);
	}

	// Start any more prefetches that need to be done.
	start_reads();

	// We successfully read another part.
	return true;
}

template<class ManifestIterator>
size_t multi_reader<ManifestIterator>::read(char* buf, size_t len)
{
	size_t tot_read = 0;
	// while we still have space to read to
	// and there's something to read from
	while (len && m_reader) {
		size_t num_read = m_reader->read(buf, len);
		if (num_read < len) {
			// SPLOG_P(LOG_DEBUG, "multi_reader::read> done reading current chunk");
			if (!next()) {
				m_reader.reset();
				SPLOG_P(LOG_DEBUG, "multi_reader::read> Done reading all the chunks in the manifest");
			}
		}
		tot_read += num_read;
		buf += num_read;
		len -= num_read;
	}
	return tot_read;
}
