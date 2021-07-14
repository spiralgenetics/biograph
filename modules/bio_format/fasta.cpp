#include "modules/bio_format/fasta.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/registry.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"
#include "modules/io/msgpack_transfer.h"


REGISTER_3(importer, fasta, readable&, bool, const std::string&);
REGISTER_3(exporter, fasta, writable&, bool, const std::string&);

void fasta_importer::import(kv_sink& sink, simple_metadata& meta)
{
	SPLOG("Importing fasta");
	std::string name = "";
	std::string seq = "";
	std::string line;
	while (m_source.readline(line, k_maxline)) {
		if (line[0] == '>') {
			if (name != "" && seq != "") {
				sink.write(name, msgpack_serialize(dna_sequence(seq)));
			}
			seq = "";
			name = line.substr(1);
		}
		else {
			seq += line;
		}
	}
	if (name != "" && seq != "") {
		sink.write(name, msgpack_serialize(dna_sequence(seq)));
	}
	SPLOG("Done importing fasta");
}

void fasta_exporter::write(const std::string & key, const std::string& value)
{
	dna_sequence seq;
	msgpack_deserialize(seq, value);

	m_sink.print(">%s\n", key.c_str());
	for (size_t i = 0; i < seq.size(); i += 80) {
		size_t sz = std::min(seq.size() - i, size_t(80));
		m_sink.print("%s\n", seq.subseq(i, sz).as_string().c_str());
	}
}

