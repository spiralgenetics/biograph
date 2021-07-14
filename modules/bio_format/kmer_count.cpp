#include "modules/bio_format/kmer_count.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"
#include "modules/io/msgpack_transfer.h"

REGISTER_3(exporter, kmer_count, writable&, bool, const std::string&);

void kmer_count_exporter::write(const std::string& key, const std::string& value)
{
	kmer_t the_kmer;
	kcount_pair kcount;

	msgpack_deserialize(the_kmer, key);
	msgpack_deserialize(kcount, value);
	dna_sequence seq(the_kmer, m_kmer_size);;
	m_sink.print("%s\t%d\t%d\n", seq.as_string().c_str(), kcount.fwd, kcount.rev);
}

