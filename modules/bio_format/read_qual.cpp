
#include "modules/bio_format/read_qual.h"
#include "modules/bio_base/seq_position.h"
#include "modules/io/utils.h"
#include "modules/io/registry.h"
#include <vector>
#include "modules/io/log.h"

REGISTER_3(exporter, read_qual, writable&, bool, const std::string&);

void read_qual_exporter::write_header()
{
	SPLOG("Exporting read quality data");
}

void read_qual_exporter::write(const std::string& key, const std::string& value)
{
	//                             A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"
	const char* base_translate = "\0\5\1\5\5\5\2\5\5\5\5\5\5\4\5\5\5\5\5\3\5\5\5\5\5\5";

	uint64_t count;
	msgpack_deserialize(count, value);
	char base = key[0];
	if (base == 'E') {
		uint64_t pos = ((unsigned char) key[1]);
		m_ends.add(pos, count);
	}
	else {
		uint64_t base_num = 5;
		if (base >= 'A' && base <= 'Z') {
			base_num = (uint64_t) base_translate[base - 'A'];
		}
		uint64_t qual = ((unsigned char) key[1]);
		uint64_t pos = ((unsigned char) key[2]);
		if (base_num >= 6) {
			throw io_exception("Base number too big");
		}
		if (qual >= 256) {
			throw io_exception("Qual number too big");
		}
		if (pos >= 256) {
			throw io_exception("Pos number too big");
		}

		m_all.add(qual, count);
		m_by_base[base_num].add(qual, count);
		m_by_pos[pos].add(qual, count);
		m_by_qual[qual].add(pos, count);
	}
}

void read_qual_exporter::write_footer()
{
	m_sink.print("{\n");
	m_sink.print("  \"overall_quality\" : ");
	print_stats(m_all, m_sink);
	m_sink.print(",\n");
	//print_stats("End position", m_ends, m_sink);
	m_sink.print("  \"by_base\" : {\n");
	for (int i = 0; i < 6; i++) {
		std::string bname = ((i == 5) ? "other" : std::string("") + ("ACGTN"[i]));
		m_sink.print("    \"%s\" : ", bname.c_str());
		print_stats(m_by_base[i], m_sink);
		m_sink.print((i == 5) ? "\n" : ",\n");
	}
	m_sink.print("  },\n");
	m_sink.print("  \"by_position\" : {\n");
	bool started = false;
	for (int i = 0; i < 256; i++) {
		if (m_by_pos[i].count() == 0) {
			continue;
		}
		if (started) { 
			m_sink.print(",\n"); 
		}
		started = true;
		m_sink.print("    %d : ", i);
		print_stats(m_by_pos[i], m_sink);
	}
	m_sink.print("\n  }\n}\n");
	/*
	m_sink.print("\nBy quality\n");
	for(int i = 0; i < 256; i++)
	{
		if (m_by_qual[i].count() == 0) continue;
		print_stats(printstring("  %d", (int) i), m_by_qual[i], m_sink);
	}
	*/
	SPLOG("Exporting read quality data complete");
}
