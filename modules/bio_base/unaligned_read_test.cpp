#include <gtest/gtest.h>
#include "modules/io/log.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/test/fastq_test_utils.h"
#include "modules/test/local_context.h"
#include "modules/test/test_utils.h"


TEST(unaligned_read, read_name)
{
	unaligned_read ur;
	std::string key;
	parse_read_name("FC81GR1ABXX:7:1101:1228:1965#TGACCAAN:1", key, ur);
	ASSERT_EQ(key, "FC81GR1ABXX:7:1101:1228:1965#TGACCAAN:1");
	ASSERT_EQ(ur.pair_number, 0);
	ASSERT_EQ(ur.name_suffix, "");
	ASSERT_EQ(build_read_name(key, ur), "FC81GR1ABXX:7:1101:1228:1965#TGACCAAN:1");
	
	key.clear();
	unaligned_read ur2;
	parse_read_name("HWI-ST1124:106:C15APACXX:1:1101:1469:2170 1:N:0:NGATGT", key, ur2);
	ASSERT_EQ(key, "HWI-ST1124:106:C15APACXX:1:1101:1469:2170 1:N:0:NGATGT");
	ASSERT_EQ(ur2.pair_number, 0);
	ASSERT_EQ(ur2.name_suffix, "");
	ASSERT_EQ(build_read_name(key, ur2), "HWI-ST1124:106:C15APACXX:1:1101:1469:2170 1:N:0:NGATGT");
}


struct pair_reads_params
{
	pair_reads_params() : is_sorted(false) {}
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(is_sorted, TF_STRICT);		
	}
	bool is_sorted;
	
	void validate() const {}
};


TEST(unaligned_read, pair)
{
	make_fastq_kv("golden/pairing.fq", make_path("e_coli.kvp"));
	local_context context(1, 500000, make_path("pair_reads_test"));

	manifest unpaired_reads;
	unpaired_reads.add(file_info(path(make_path("e_coli.kvp")), 2690, 10), 0);

	manifest paired_reads = context.map_reduce("identity", "", "pair", "pair", "", unpaired_reads);
	
	SPLOG("%lu unpaired reads in.  %lu paired reads out.", unpaired_reads.get_num_records(), paired_reads.get_num_records());
	// There are two unpaired read in the fastq, so we have to subtract them out.
	ASSERT_EQ(unpaired_reads.get_num_records(), 2 * (paired_reads.get_num_records() - 2));
	
	manifest_reader reads_reader(paired_reads);
	read_id key;
	unaligned_reads value;
	unsigned paired_count = 0;
	unsigned unpaired_count = 0;
	while(reads_reader.read_msgpack(key, value)) {
		ASSERT_TRUE((value.size() == 1) || (value.size() == 2));
		if (value.size() == 2) { paired_count++; }
		if (value.size() == 1) { unpaired_count++; }
	}
	ASSERT_EQ(paired_count, 5);
	ASSERT_EQ(unpaired_count, 2);
}


TEST(unaligned_read, pair_no_suffix)
{
	make_fastq_kv("golden/pairing_no_suffix.fq", make_path("e_coli.kvp"));
	local_context context(1, 500000, make_path("no_suffix_pair_reads_test"));

	manifest unpaired_reads;
	unpaired_reads.add(file_info(path(make_path("e_coli.kvp")), 2690, 10), 0);

	manifest paired_reads = context.map_reduce("identity", "", "lexical", "pair", "", unpaired_reads);
	
	SPLOG("%lu unpaired reads in.  %lu paired reads out.", unpaired_reads.get_num_records(), paired_reads.get_num_records());
	// There are two unpaired read in the fastq, so we have to subtract them out.
	ASSERT_EQ(unpaired_reads.get_num_records(), 2 * (paired_reads.get_num_records() - 2));
	
	manifest_reader reads_reader(paired_reads);
	read_id key;
	unaligned_reads value;
	unsigned paired_count = 0;
	unsigned unpaired_count = 0;
	while(reads_reader.read_msgpack(key, value)) {
		ASSERT_TRUE((value.size() == 1) || (value.size() == 2));
		if (value.size() == 2) { paired_count++; }
		if (value.size() == 1) { unpaired_count++; }
	}
	ASSERT_EQ(paired_count, 5);
	ASSERT_EQ(unpaired_count, 2);
}
