#include "gtest/gtest.h"
#include "modules/test/test_utils.h"
#include "modules/test/fastq_test_utils.h" 
#include "modules/bio_base/unaligned_read.h" 
#include "modules/bio_format/fastq.h"
#include "modules/io/file_io.h"
#include "modules/io/keyvalue.h"

static const int line_size = 1<<16;

TEST(fastq, import)
{
	// Import FASTQ into MSGPACK-serialized kvp file.
	ASSERT_NO_THROW(make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp")));

	// Verify that line 1,2 and 4 made it into an unaligned_read structure.
	file_reader original_fastq("golden/e_coli_10000snp.fq");
	file_reader spiral_fastq(make_path("e_coli_10000.kvp"));
	kv_reader spiral_fastq_as_kv(spiral_fastq);

	read_id seq_id;
	unaligned_reads paired_reads;

	while (spiral_fastq_as_kv.read_msgpack (seq_id, paired_reads)) {
		for (unaligned_read one_read : paired_reads) {

			std::string sequence_identifier;
			std::string dna_sequence;
			std::string plus_comments;
			std::string quality_sequence;

			ASSERT_TRUE(original_fastq.readline(sequence_identifier, line_size));
			ASSERT_TRUE(original_fastq.readline(dna_sequence, line_size));
			ASSERT_TRUE(original_fastq.readline(plus_comments, line_size));
			ASSERT_TRUE(original_fastq.readline(quality_sequence, line_size));

			// we don't store the '@' character
			ASSERT_EQ(one_read.original_sequence_id.compare(sequence_identifier.substr(1)), 0); 
			ASSERT_EQ(one_read.sequence.compare(dna_sequence), 0);
			ASSERT_EQ(one_read.quality.compare(quality_sequence), 0);
		}
	}
}

TEST(fastq, export_)
{
	// Import FASTQ into MSGPACK-serialized kvp file.
	ASSERT_NO_THROW( make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp")) );

	// Export the FASTQ file
	file_reader in_kvp(make_path("e_coli_10000.kvp"));
	file_writer out_fastq(make_path("e_coli_10000.fastq"));
	kv_reader kin(in_kvp);
	fastq_exporter a_fastq_exporter(out_fastq);

	ASSERT_NO_THROW(a_fastq_exporter.export_from(kin));

	in_kvp.close();
	out_fastq.close();

	// Compare the exported FASTQ file with the original
	std::string original_line, spiral_line;
	file_reader original_fastq("golden/e_coli_10000snp.fq");
	file_reader spiral_fastq(make_path("e_coli_10000.fastq"));
	bool continue_reading = true;
	int num_blobs = 0;
	while (continue_reading) {
		// 1st line is sequence identifier
		continue_reading = original_fastq.readline( original_line, line_size);
		if (!continue_reading) {
			break;
		}
		ASSERT_TRUE( spiral_fastq.readline( spiral_line, line_size) );
		ASSERT_EQ( spiral_line.compare( original_line ), 0);

		// 2nd line is dna sequence
		ASSERT_TRUE( spiral_fastq.readline( spiral_line, line_size) );
		ASSERT_TRUE( original_fastq.readline( original_line, line_size) );
		ASSERT_EQ( spiral_line.compare( original_line ), 0);

		// 3rd line is + comments
		ASSERT_TRUE( spiral_fastq.readline( spiral_line, line_size) );
		ASSERT_TRUE( original_fastq.readline( original_line, line_size) );
		// They don't need to match since we don't store them in Spiral

		// 4th line is quality
		ASSERT_TRUE( spiral_fastq.readline( spiral_line, line_size) );
		ASSERT_TRUE( original_fastq.readline( original_line, line_size) );
		ASSERT_EQ( spiral_line.compare( original_line ), 0);

		num_blobs++;
	}
	ASSERT_EQ( num_blobs, 10000 );
}
