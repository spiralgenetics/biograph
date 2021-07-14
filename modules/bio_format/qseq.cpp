#include "modules/bio_format/qseq.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"
#include "modules/io/registry.h"

#include <vector>

static const int k_maxline = 1000;

REGISTER_3(importer, qseq, readable&, bool, const std::string&);

void qseq_importer::import(kv_sink& sink, simple_metadata& meta)
{
	SPLOG("Importing qseq");
	std::string line;
	int linenum = 0;
	size_t bases = 0;
	while (m_source.readline(line, k_maxline)) {
		read_id id;
		unaligned_read read;
		unaligned_reads reads;
		linenum++;
		std::vector<std::string> parts;
		size_t start = 0;
		for (size_t i = 0; i < line.size(); i++) {
			if (line[i] == '\t') {
				parts.push_back(line.substr(start, i - start));
				start = i + 1;
			}
		}
		if (parts.size() < 10) {
			throw io_exception(printstring(
				"Line %d: Not enough tab delimited columns", linenum));
		}

		id.pair_name = parts[0] + "_" 
				 + parts[1] + ":" 
				 + parts[2] + ":" 
				 + parts[3] + ":" 
				 + parts[4] + ":" 
				 + parts[5] + ":"
				 + parts[6];
		read.pair_number = atoi(parts[7].c_str());
		std::string seq = parts[8];
		std::string qual = parts[9];
		if (seq.size() != qual.size()) {
			throw io_exception(printstring(
				"Line %d: Size of quality and nuclitide data doesn't match", linenum));
		}

		for (size_t i = 0; i < seq.size(); i++) {
			if (seq[i] == '.' || qual[i] == 'B') {
				seq[i] = 'N';
			}
			if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T' && seq[i] != 'N') {
				throw io_exception(printstring(
					"Line %d: Invalid char in sequence", linenum));
			}
			if (qual[i] < 64) {
				throw io_exception(printstring(
					"Line %d: Invalid char in quality data", linenum));
			}
			qual[i] -= 31;
		}

		bool found_first = false;
		size_t first = 0;
		size_t last = 0;
		for (size_t i = 0; i < seq.size(); i++) {
			if (seq[i] != 'N') {
				if (found_first == false) {
					first = i;
					found_first = true;
				}
				last = i;
			}
		}
		
		size_t len = last - first + 1;

		if (found_first && len >= 30) {
			read.sequence = seq.substr(first, len);
			read.quality = qual.substr(first, len);
			bases += read.sequence.size();
			reads.push_back(read);
			sink.write_msgpack(id, reads);
		}
	}

	SPLOG("Done importing qseq");
	meta.set_simple("sample_bases", bases);
}
