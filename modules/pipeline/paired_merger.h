#pragma once

#include "modules/bio_base/unaligned_read.h"
#include "modules/io/keyvalue.h"

class paired_merger : public kv_sink
{
public:
	paired_merger(kv_sink& out, kv_source& merge)
		: m_out(out)
		, m_merge(merge)
	{}

	void write(const std::string& key, const std::string& value) override
	{
		read_id key1, key2;
		unaligned_reads reads1, reads2;
		// Convert incoming data from strings to reads
		msgpack_deserialize(key1, key);
		msgpack_deserialize(reads1, value);
		// Read data from existing file
		bool more = m_merge.read_msgpack(key2, reads2);
		// Make sure old file has data
		if (!more) {
			//throw io_exception("During pairing first file ran out of records at key '" + key + "'");
			unaligned_reads both;
			both.push_back(reads1[0]);
			m_out.write_msgpack(key1, both);
		} else {
			// Check that both are single reads, and that they are 1/2
			// Maybe we shouldn't care about the order (ask Becky) [we do care about the order (~Adam)]
			if (reads1.size() != 1) {
				throw io_exception(printstring(
					"In 2nd paired file, unexpected reads size of %lu, should have been 1 for key: %s",
					reads1.size(), key1.pair_name.c_str()).c_str()
				);
			}
			if (reads2.size() != 1) {
				throw io_exception(printstring(
					"In 1st paired file, unexpected reads size of %lu, should have been 1 for key: %s",
					reads1.size(), key2.pair_name.c_str()).c_str()
				);
			}
	
			// Combine the reads
			unaligned_reads both;
			both.push_back(reads1[0]);
			both.push_back(reads2[0]);
			// Write to the downstream output.
			m_out.write_msgpack(key1, both);
		}
	}

	void close() override
	{
		std::string key2;
		std::string value2;
		bool more = m_merge.read(key2, value2);
		if (more) {
			throw io_exception("During pairing second file ran out of records at key '" + key2 + "'");
		}
		m_out.close();
	}

private:
	kv_sink& m_out;
	kv_source& m_merge;
};
