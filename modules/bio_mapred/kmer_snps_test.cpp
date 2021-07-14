
#include "modules/bio_mapred/kmer_snps.h"
#include "modules/mapred/kv_hold.h"
#include "modules/io/log.h"
#include <gtest/gtest.h>

std::string kmer_serial = 
"{\"kmer_size\":30,\"lookup\":{\"all_metadata\":{\"m_data\":{}},\"meta\":\"\",\"num_records\":0,\"partitions\":[{\"files\":[{\"file\":\"yeast_tmp/5ef1ec91-b195-49f5-86e6-1cf82950c2ba\",\"first_key\":\"\",\"last_key\":\"\",\"num_records\":0,\"size\":4194312}],\"num_records\":0,\"size\":4194312}],\"size\":4194312,\"sort\":\"\",\"tags\":{}},\"size\":11749598,\"table\":{\"all_metadata\":{\"m_data\":{}},\"meta\":\"\",\"num_records\":0,\"partitions\":[{\"files\":[{\"file\":\"yeast_tmp/7c44d20a-45f0-4d3a-a6ca-d13977789a93\",\"first_key\":\"\",\"last_key\":\"\",\"num_records\":0,\"size\":58747990}],\"num_records\":0,\"size\":58747990}],\"size\":58747990,\"sort\":\"\",\"tags\":{}}}";

/*
"{\"kmer_size\":30,\"lookup\":{\"all_metadata\":{\"m_data\":{}},\"meta\":\"\",\"num_records\":0,\"partitions\":[{\"files\":[{\"file\":\"yeast_tmp/c3ab5020-5f25-4535-a190-82ef43967480\",\"first_key\":\"\",\"last_key\":\"\",\"num_records\":0,\"size\":4194312}],\"num_records\":0,\"size\":4194312}],\"size\":4194312,\"sort\":\"\",\"tags\":{}},\"size\":12033652,\"table\":{\"all_metadata\":{\"m_data\":{}},\"meta\":\"\",\"num_records\":0,\"partitions\":[{\"files\":[{\"file\":\"yeast_tmp/6569dc0b-e9a2-4405-8460-381436213c8d\",\"first_key\":\"\",\"last_key\":\"\",\"num_records\":0,\"size\":60168260}],\"num_records\":0,\"size\":60168260}],\"size\":60168260,\"sort\":\"\",\"tags\":{}}}";
*/
void kmer_write(kv_sink& s, const char* kmer) 
{
	kmer_t k = dna_sequence(kmer).as_kmer();
	s.write(msgpack_serialize(k), "X");
}

TEST(kmers_snps, basic) 
{
	SPLOG("Hey!");
	kv_hold h("");
	kmer_write(h, "AAAAGGG");
	kmer_write(h, "ACCTGAG");
	kmer_write(h, "CCAGAGA");
	kmer_write(h, "CCCTTTA");
	kmer_write(h, "CCTGAGA");
	kmer_set ks(h, 5, 7);
	std::map<kmer_t, kmer_t> out_map;
	kmer_find_snps(ks, out_map, 4000000000l, 16);
	SPLOG("%d", int(ks.size()));
}

/*
TEST(kmers_snps, magic) 
{
	SPLOG("Hey!");
	kmer_set ks(kmer_serial);
	std::map<kmer_t, kmer_t> out_map;
	kmer_find_snps(ks, out_map, 4000000000l, 16);
	SPLOG("%d", int(ks.size()));
}
*/

