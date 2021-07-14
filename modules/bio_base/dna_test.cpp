#include "modules/bio_base/kmer.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/dna_base_set.h"
#include "modules/bio_base/dna_multiseq.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/log.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <chrono>
#include <boost/format.hpp>

dna_sequence make_random_sequence(unsigned int size)
{
	std::string dna_string;
	dna_string.reserve(size);
	for(unsigned int i = 0; i < size; i++)
	{
		dna_string.push_back("ACGT"[random() % 4]);
	}
	
	return dna_sequence(dna_string);
}

TEST(dna, base)
{
	dna_base b('C');
	int x = int(b);
	ASSERT_EQ(x, 1);
}

TEST(dna, all_bases)
{
  dna_base_array<dna_base> base_array;
  dna_base_array<int> int_array;
  for (dna_base b : dna_bases()) {
    base_array[b] = b;
    int_array[b] = 10 + int(b);
  }
  EXPECT_THAT(base_array, testing::ElementsAre(dna_base(0), dna_base(1),
                                              dna_base(2), dna_base(3)));
  EXPECT_THAT(int_array, testing::ElementsAre(10, 11, 12, 13));
}

TEST(dna, encoding)
{
	size_t size;
	std::string pre;
	for(size_t i = 0; i < 1000; i++)
	{
		size = random() % 500;
		pre.resize(size);
		for(size_t i = 0; i < size; i++)
			pre[i] = "ACGT"[random() % 4];
		dna_sequence seq(pre);
		ASSERT_EQ(seq.size(), size);
		std::string post = seq.as_string();
		ASSERT_EQ(pre, post);
		std::string bin = seq.as_packed();
		ASSERT_EQ(bin.size(), size/4 + 1);
		dna_sequence seq2(bin, true);
		std::string post2 = seq2.as_string();
		ASSERT_EQ(pre, post2);
	}
	size = 9;
	pre.resize(size);
	for(size_t i = 0; i < size; i++)
		pre[i] = "ACGT"[random() % 4];
	printf("%s\n", pre.c_str());
	dna_sequence s(pre);
	std::string bin = s.as_packed();
	for(size_t i = 0; i < bin.size(); i++)
	{
		unsigned char c = bin[i];
		for(int j = 7; j >= 0; j--)
			printf("%d", (c >> j) & 1);
		printf(" ");
	}
	printf("\n");
}

TEST(dna, complement)
{
	dna_sequence a("AATGTAGCCTAG");
	std::string ca("CTAGGCTACATT");
	ASSERT_EQ(a.rev_comp().as_string(), ca);
}

TEST(dna, base_set)
{
	dna_base_set a('A');
	dna_base_set c('C');
	dna_base_set ac = a | c;
	ASSERT_EQ(ac.as_list(), "A,C");
	ASSERT_EQ(ac.as_code(), 'M');
}

TEST(dna, mutate)
{
	size_t size;
	std::string str;
	for(size_t i = 0; i < 10000; i++)
	{
		size = random() % 50 + 1;
		str.resize(size);
		for(size_t i = 0; i < size; i++)
			str[i] = "ACGT"[random() % 4];
		dna_sequence seq1(str);
		dna_sequence seq2(seq1);
		dna_sequence seq(size);
		for(size_t i = 0; i < size; i++)
			seq[i] = seq2[i];
		// Set one
		char nval = "ACGT"[random() % 4];
		size_t loc = random() % size;
		seq[loc] = dna_base(nval);
		str[loc] = nval;
		// Move one
		size_t loc1 = random() % size;
		size_t loc2 = random() % size;
		seq[loc1] = seq[loc2];
		str[loc1] = str[loc2];
		std::string seqstr = seq.as_string();
		ASSERT_EQ(str, seqstr);
	}
}

TEST(dna, multiseq)
{
	// ...................**......<..............>...........
	dna_sequence a("ATAGCCAAGATCCAAGGACCTAGTCATCGTACACAAGCCA");
	dna_sequence b("TAGCCGGGATCCATAGGACCTAGTCATCTACACAAGCCAT");

	dna_multiseq ms(a,b);

	std::string ao = ms.get_string(0);
	std::string bo = ms.get_string(1);
	printf("%s\n%s\n", ao.c_str(), bo.c_str());
	
	ASSERT_EQ(ao, "ATAGCCAAGATCCA.AGGACCTAGTCATCGTACACAAGCCA.");
	ASSERT_EQ(bo, ".TAGCCGGGATCCATAGGACCTAGTCATC.TACACAAGCCAT");
}

TEST(dna, kmer)
{

	//            1234567890123456789
	dna_sequence s("TCCTAGTAGACATGCCATG");
	size_t k = s.as_kmer();
	dna_sequence s2 = dna_sequence(k, 19);
	printf("%s\n", s2.as_string().c_str());
	ASSERT_EQ(s.as_string(), s2.as_string());
	int b = int(dna_base('C'));
	rotate_left(k, 19, b);
	printf("%s\n", dna_sequence(k, 19).as_string().c_str());
	ASSERT_EQ(dna_sequence(k, 19).as_string(), "CCTAGTAGACATGCCATGC");
	ASSERT_EQ(b, (int) dna_base('T'));
	rotate_right(k, 19, b);
	printf("%s\n", dna_sequence(k, 19).as_string().c_str());
	ASSERT_EQ(dna_sequence(k, 19).as_string(), "TCCTAGTAGACATGCCATG");
	ASSERT_EQ(b, (int) dna_base('C'));

	printf("0x%016zx\n", dna_sequence("A").as_kmer());
	printf("0x%016zx\n", dna_sequence("AA").as_kmer());
	printf("0x%016zx\n", dna_sequence("AAAA").as_kmer());
	printf("0x%016zx\n", dna_sequence("C").as_kmer());
	printf("0x%016zx\n", dna_sequence("CC").as_kmer());
	printf("0x%016zx\n", dna_sequence("CCCC").as_kmer());
	printf("0x%016zx\n", dna_sequence("G").as_kmer());
	printf("0x%016zx\n", dna_sequence("GG").as_kmer());
	printf("0x%016zx\n", dna_sequence("GGGG").as_kmer());
	printf("0x%016zx\n", dna_sequence("T").as_kmer());
	printf("0x%016zx\n", dna_sequence("TT").as_kmer());
	printf("0x%016zx\n", dna_sequence("TTTT").as_kmer());
}

TEST(dna, kmer2)
{
	for(size_t ks = 1; ks < 19; ks++)
	{
		for(size_t tc = 0; tc < 20; tc++)
		{
			dna_sequence s;
			for(size_t i = 0; i < ks; i++)
				s.push_back(dna_base(rand()%4));
			kmer_t k = s.as_kmer();
			dna_sequence s2 = dna_sequence(k, ks);
			//SPLOG("orig = %s, kmer = %d, new = %s\n", s.as_string().c_str(), (int) k, s2.as_string().c_str());
			ASSERT_EQ(s.as_string(), s2.as_string());
		}
	}	
}

TEST(dna, kmerize)
{
	const size_t kmer_size = 9;
	dna_sequence seq("ACGTACGTACGTACGTACGT");
	SPLOG("%s", seq.as_string().c_str());

	std::vector<std::string> kmer_strs;
	std::vector<kmer_t> kmers;

	kmer_t kmer = make_kmer(seq.begin(), kmer_size);
	kmers.push_back(kmer);
	SPLOG("0x%016zx", kmer);
	kmer_strs.push_back(dna_sequence(kmer, kmer_size).as_string());
	for(size_t offset = 1; offset <= seq.size() - kmer_size; offset++) {
		int base = int(seq[kmer_size + offset - 1]);
		rotate_left(kmer, kmer_size, base);
		kmer_strs.push_back(dna_sequence(kmer, kmer_size).as_string());
		kmers.push_back(kmer);
		SPLOG("0x%016zx", kmer);
	}

	ASSERT_EQ(seq.size() - kmer_size + 1, kmer_strs.size());
	ASSERT_EQ("ACGTACGTA", kmer_strs[0]); 
	ASSERT_EQ("CGTACGTAC", kmer_strs[1]);
	ASSERT_EQ("GTACGTACG", kmer_strs[2]);
	ASSERT_EQ("TACGTACGT", kmer_strs[3]);
	ASSERT_EQ("ACGTACGTA", kmer_strs[4]);
	ASSERT_EQ("CGTACGTAC", kmer_strs[5]);
	ASSERT_EQ("GTACGTACG", kmer_strs[6]);
	ASSERT_EQ("TACGTACGT", kmer_strs[7]);
	ASSERT_EQ("ACGTACGTA", kmer_strs[8]);
	ASSERT_EQ("CGTACGTAC", kmer_strs[9]);
	ASSERT_EQ("GTACGTACG", kmer_strs[10]);
	ASSERT_EQ("TACGTACGT", kmer_strs[11]);

	ASSERT_EQ(0x0000000000006c6c, kmers[0]);
	ASSERT_EQ(0x000000000001b1b1, kmers[1]);
	ASSERT_EQ(0x000000000002c6c6, kmers[2]);
	ASSERT_EQ(0x0000000000031b1b, kmers[3]);
	ASSERT_EQ(0x0000000000006c6c, kmers[4]);
	ASSERT_EQ(0x000000000001b1b1, kmers[5]);
	ASSERT_EQ(0x000000000002c6c6, kmers[6]);
	ASSERT_EQ(0x0000000000031b1b, kmers[7]);
	ASSERT_EQ(0x0000000000006c6c, kmers[8]);
	ASSERT_EQ(0x000000000001b1b1, kmers[9]);
	ASSERT_EQ(0x000000000002c6c6, kmers[10]);
	ASSERT_EQ(0x0000000000031b1b, kmers[11]);
}

TEST(dna, reverse)
{
	dna_sequence s("TCCTAGTAGACATGCCATG");
	std::string tot;
	for(auto i = s.rcbegin(); i != s.rcend(); i++)
	{
		char base = char(*i);
		tot += base;
	}
	printf("%s\n", tot.c_str());
	ASSERT_EQ(tot, "CATGGCATGTCTACTAGGA");
}

TEST(dna, byte_reverse)
{
	ASSERT_EQ(byte_rev_comp_bases(0), 0xff);
	ASSERT_EQ(byte_rev_comp_bases(0xff), 0);
	ASSERT_EQ(byte_rev_comp_bases(0xd5), 0xa8);
	ASSERT_EQ(byte_rev_comp_bases(0x12), 0x7b);
	ASSERT_EQ(byte_rev_comp_bases(0xc3), 0x3c);
	ASSERT_EQ(byte_rev_comp_bases(0x44), 0xee);
	ASSERT_EQ(byte_rev_comp_bases(0x3d), 0x83);
	ASSERT_EQ(byte_rev_comp_bases(0xe8), 0xd4);
}

TEST(dna, basic_equality)
{
	dna_sequence seq1;
	dna_sequence seq2;
	
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq1.begin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.begin(), seq1.size()));
	seq1 = dna_sequence("A");
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq1.begin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.begin(), seq2.size()));
	ASSERT_TRUE(subseq_equal(seq2.begin(), seq1.begin(), seq2.size()));
	seq2 = dna_sequence("A");
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.begin(), seq1.size()));
	seq2 = dna_sequence("T");
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.rcbegin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq1.rcbegin(), seq2.begin(), seq1.size()));
	seq2 = dna_sequence("G");
	ASSERT_FALSE(subseq_equal(seq1.begin(), seq2.begin(), seq1.size()));
	seq1 = dna_sequence("G");
	seq2 = dna_sequence("GG");
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.begin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq2.begin(), seq1.begin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.begin() + 1, seq1.size()));
	ASSERT_TRUE(subseq_equal(seq2.begin() + 1, seq1.begin(), seq1.size()));
	seq1 = dna_sequence("TA");
	seq2 = dna_sequence("TCTA");
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.begin() + 2, seq1.size()));
	ASSERT_FALSE(subseq_equal(seq1.begin(), seq2.begin() + 1, seq1.size()));
	seq1 = dna_sequence("TAGACCTGCCGGATATAA");
	seq2 = dna_sequence("CCGTATGATAGCCGTAGG");
	ASSERT_TRUE(subseq_equal(seq1.begin() + 7, seq2.begin() + 10, 4));
	ASSERT_FALSE(subseq_equal(seq1.begin(), seq2.begin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq1.rcbegin() + 6, seq2.begin(), 3));
	seq2 = dna_sequence("TTATATCCGGCAGGTCTA");
	ASSERT_TRUE(subseq_equal(seq1.rcbegin(), seq2.begin(), seq1.size()));
	ASSERT_TRUE(subseq_equal(seq1.begin(), seq2.rcbegin(), seq1.size()));
	ASSERT_FALSE(subseq_equal(seq1.begin(), seq2.begin(), seq1.size()));
}

TEST(dna, subseq_compare)
{
	unsigned int seed = static_cast<unsigned int>(::random());
	SPLOG("Seed = %u", seed);
	::srandom(seed);
	for (unsigned int i = 0; i < 10000; i++)
	{
		unsigned int seq_size = ::random() % 500 + 1;
		unsigned int subseq_size = ::random() % seq_size;
		unsigned int offset = ::random() % (seq_size - subseq_size);
		bool complement = ::random() % 2;
		dna_sequence seq(make_random_sequence(seq_size));
		if (complement)
		{
			dna_sequence subseq(seq.subseq(offset, subseq_size).rev_comp());
			for (unsigned int j = subseq.size(); j > 0; j--)
			{
				ASSERT_TRUE(subseq_equal(seq.begin() + offset, subseq.rcbegin(), j));
			}
			
			for (unsigned int j = subseq.size(); j > 0; j--)
			{
				ASSERT_TRUE(subseq_equal(seq.begin() + offset + j, subseq.rcbegin() + j, subseq.size() - j));
			}
			
			if (subseq.size() == 0) continue;
			unsigned int random_change = ::random() % subseq.size();
			dna_base random_base = subseq[random_change];
			subseq[random_change] = random_base.complement();
			ASSERT_FALSE(subseq_equal(seq.begin() + offset, subseq.rcbegin(), subseq.size()));
		}
		else
		{
			dna_sequence subseq(seq.subseq(offset, subseq_size));
			for (unsigned int j = subseq.size(); j > 0; j--)
			{
				ASSERT_TRUE(subseq_equal(seq.begin() + offset, subseq.begin(), j));
			}
			
			for (unsigned int j = subseq.size(); j > 0; j--)
			{
				ASSERT_TRUE(subseq_equal(seq.begin() + offset + j, subseq.begin() + j, subseq.size() - j));
			}
			
			if (subseq.size() == 0) continue;
			unsigned int random_change = ::random() % subseq.size();
			dna_base random_base = subseq[random_change];
			subseq[random_change] = random_base.complement();
			ASSERT_FALSE(subseq_equal(seq.begin() + offset, subseq.begin(), subseq.size()));
		}
	}
}

TEST(dna, basic_lessthan)
{
	dna_sequence seq1;
	dna_sequence seq2;
	
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq2.begin(), seq1.begin(), seq2.size(), seq1.size()));
	seq1 = dna_sequence("A");
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq1.begin(), seq1.size(), seq1.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_TRUE(subseq_lessthan(seq2.begin(), seq1.begin(), seq2.size(), seq1.size()));
	seq2 = dna_sequence("A");
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	seq2 = dna_sequence("T");
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.rcbegin(), seq1.size(), seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.rcbegin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_TRUE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq2.begin(), seq1.begin(), seq2.size(), seq1.size()));
	seq2 = dna_sequence("G");
	ASSERT_TRUE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	seq1 = dna_sequence("G");
	seq2 = dna_sequence("GG");
	ASSERT_TRUE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq2.begin(), seq1.begin(), seq2.size(), seq1.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin() + 1, seq1.size(), seq2.size() - 1));
	ASSERT_FALSE(subseq_lessthan(seq2.begin() + 1, seq1.begin(), seq2.size() - 1, seq1.size()));
	seq1 = dna_sequence("TAT");
	seq2 = dna_sequence("TCTAA");
	ASSERT_TRUE(subseq_lessthan(seq2.begin() + 2, seq1.begin(), 3, 3));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin() + 2, 3, 3));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin() + 2, 2, 2));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin() + 1, seq1.size(), seq2.size() - 1));
	seq1 = dna_sequence("TAGACCTGCCGGATATAA");
	seq2 = dna_sequence("CCGTATGATAGCCGTAGG");
	ASSERT_FALSE(subseq_lessthan(seq1.begin() + 7, seq2.begin() + 10, 4, 4));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.rcbegin(), seq2.rcbegin(), seq1.size(), seq2.size()));
	ASSERT_TRUE(subseq_lessthan(seq2.begin(), seq1.begin(), seq2.size(), seq1.size()));
	ASSERT_TRUE(subseq_lessthan(seq2.rcbegin(), seq1.rcbegin(), seq2.size(), seq1.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.rcbegin() + 6, seq2.begin(), 3, 3));
	ASSERT_TRUE(subseq_lessthan(seq1.rcbegin() + 6, seq2.begin(), 3, 4));
	seq2 = dna_sequence("TTATATCCGGCAGGTCTA");
	ASSERT_FALSE(subseq_lessthan(seq1.rcbegin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.begin(), seq2.rcbegin(), seq1.size(), seq2.size()));
	ASSERT_TRUE(subseq_lessthan(seq1.begin(), seq2.begin(), seq1.size(), seq2.size()));
	ASSERT_TRUE(subseq_lessthan(seq1.rcbegin(), seq2.begin(), seq1.size() - 1, seq2.size()));
	ASSERT_FALSE(subseq_lessthan(seq1.rcbegin(), seq2.begin(), seq1.size(), seq2.size() - 1));
}

struct dna_holder
{
	dna_sequence m_seq;
	ptrdiff_t m_offset;
	size_t m_length; // Comparison length, not sequence length.
	bool m_rev_comp;
	
	dna_holder();
	explicit dna_holder(size_t a_constant_size);
};

dna_holder::dna_holder()
{
	unsigned int seq_size = ::random() % 20 + 1;
	m_offset = ::random() % seq_size;
	m_length = ::random() % (seq_size - m_offset);
	m_rev_comp = ::random() % 2;
	m_seq = make_random_sequence(seq_size);
}

dna_holder::dna_holder(size_t a_constant_size)
{
	m_offset = ::random() % a_constant_size;
	m_length = ::random() % (a_constant_size - m_offset);
	m_rev_comp = ::random() % 2;
	m_seq = make_random_sequence(a_constant_size);
}

inline bool operator<(const dna_holder& lhs, const dna_holder& rhs)
{
	dna_const_iterator lhs_it = lhs.m_rev_comp ? lhs.m_seq.rcbegin() : lhs.m_seq.begin();
	dna_const_iterator rhs_it = rhs.m_rev_comp ? rhs.m_seq.rcbegin() : rhs.m_seq.begin();
	return dna_sequence(lhs_it + lhs.m_offset, lhs_it + lhs.m_offset + lhs.m_length) <
		dna_sequence(rhs_it + rhs.m_offset, rhs_it + rhs.m_offset + rhs.m_length);
	return subseq_lessthan(lhs_it + lhs.m_offset, rhs_it + rhs.m_offset, lhs.m_length, rhs.m_length);
}

TEST(dna, sort)
{
	unsigned int seed = static_cast<unsigned int>(::random());
	SPLOG("Seed = %u", seed);
	::srandom(seed);
	
	std::vector<dna_holder>	seqs(100000);
	
	std::sort(seqs.begin(), seqs.end());
	
	dna_sequence previous_seq;
	dna_sequence current_seq;
	for (const auto& seq : seqs)
	{
		if (seq.m_rev_comp) current_seq = seq.m_seq.rev_comp().subseq(seq.m_offset, seq.m_length);
		else current_seq = seq.m_seq.subseq(seq.m_offset, seq.m_length);
		ASSERT_TRUE(previous_seq.as_string() <= current_seq.as_string());
		previous_seq = current_seq;
	}
}

TEST(dna, sort_same_size)
{
	unsigned int seed = static_cast<unsigned int>(::random());
	SPLOG("Seed = %u", seed);
	::srandom(1637635076);
	
	std::vector<dna_holder>	seqs;
	const unsigned int  k_size = 100000;
	seqs.reserve(k_size);
	for (unsigned int i = 0; i < k_size; i++)
	{
		seqs.emplace_back(dna_holder(100));
	}
	
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
	std::sort(seqs.begin(), seqs.end());
    end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end-start;
	SPLOG("Sorted %u dna_sequences in %f seconds", k_size, elapsed.count());
	
	dna_sequence previous_seq;
	dna_sequence current_seq;
	for (const auto& seq : seqs)
	{
		if (seq.m_rev_comp) current_seq = seq.m_seq.rev_comp().subseq(seq.m_offset, seq.m_length);
		else current_seq = seq.m_seq.subseq(seq.m_offset, seq.m_length);
		ASSERT_TRUE(previous_seq.as_string() <= current_seq.as_string());
		previous_seq = current_seq;
	}
}

TEST(dna, reserve)
{
	dna_sequence seq;
	seq.reserve(1000);
	seq.push_back(dna_base('A'));
	dna_iterator small_begin = seq.begin();
	for (unsigned int i = 0; i < 999; i++)
	{
		seq.push_back(dna_base('T'));
	}
	dna_iterator large_begin = seq.begin();
	ASSERT_EQ(small_begin.get_original_data(), large_begin.get_original_data());
	seq.push_back(dna_base('C'));
	ASSERT_NE(small_begin.get_original_data(), seq.begin().get_original_data());
}

TEST(dna, iupac_rc)
{
	std::string test_string{"YWVTSRNMKHGDCBAABCDHKMNRSTVWY"};
	reverse_complement_iupac_string(test_string);
	ASSERT_EQ(test_string, "RWBASYNKMDHGVTTVGHCDMKNYSABWR");
	
	test_string = "AATGTAGCCTAG";
	reverse_complement_iupac_string(test_string);
	ASSERT_EQ(test_string, "CTAGGCTACATT");
	
	test_string = "NNNNNNNNNNNN";
	reverse_complement_iupac_string(test_string);
	ASSERT_EQ(test_string, "NNNNNNNNNNNN");
}

class dna_compare_test : public testing::Test {
 public:
  static constexpr size_t k_max_compare_size = 28 * 3;
  enum action_t { NONE, CHANGE_TO_BIG, CHANGE_TO_SMALL, TRUNCATE };

  dna_slice make_slice(unsigned initial_offset, action_t action,
                       unsigned action_pos, bool rc) {
    dna_sequence prefix_seq = dna_sequence("AAAA").subseq(0, initial_offset);

    dna_sequence seq;
    for (size_t i = 0; i < k_max_compare_size; ++i) {
      seq.push_back(dna_base((i & 1) ? 'C' : 'G'));
    }
    switch (action) {
      case NONE:
        break;
      case CHANGE_TO_BIG:
        seq[action_pos] = dna_base('T');
        break;
      case CHANGE_TO_SMALL:
        seq[action_pos] = dna_base('A');
        break;
      case TRUNCATE:
        seq = seq.subseq(0, action_pos);
        break;
    }
    if (rc) {
      m_seqs.push_back(prefix_seq + seq);
      return dna_slice(m_seqs.back().begin() + prefix_seq.size(), seq.size());
    } else {
      m_seqs.push_back((prefix_seq + seq).rev_comp());
      return dna_slice(m_seqs.back().rcbegin() + prefix_seq.size(), seq.size());
    }
  }

  void check_shared(const dna_slice& slice1, const dna_slice& slice2) {
    unsigned expected_shared = 0;
    while (expected_shared < slice1.size() && expected_shared < slice2.size() &&
           slice1[expected_shared] == slice2[expected_shared]) {
      expected_shared++;
    }

    EXPECT_EQ(expected_shared, slice1.shared_prefix_length(slice2))
        << "\nSlice1: " << slice1 << "\nSlice2: " << slice2;
    EXPECT_EQ(expected_shared, slice2.shared_prefix_length(slice1))
        << "\nSlice1: " << slice1 << "\nSlice2: " << slice2;
  }

 private:
  std::list<dna_sequence> m_seqs;
};

TEST_F(dna_compare_test, equal) {
  for (unsigned offset1 = 0; offset1 < 4; ++offset1) {
    for (unsigned offset2 = 0; offset2 < 4; ++offset2) {
      for (bool rc1 : {false, true}) {
        for (bool rc2 : {false, true}) {
          dna_slice slice1 = make_slice(offset1, NONE, 0, rc1);
          dna_slice slice2 = make_slice(offset2, NONE, 0, rc2);
          EXPECT_EQ(dna_compare_result::EQUAL, slice1.compare_to(slice2))
              << "Offset1: " << offset1 << " Offset2: " << offset2
              << " Rc1: " << rc1 << " Rc2: " << rc2 << " Slice1:\n"
              << slice1.as_string() << " Slice2:\n"
              << slice2.as_string();
            check_shared(slice1, slice2);
        }
      }
    }
  }
}

TEST_F(dna_compare_test, less_than) {
  for (unsigned offset1 = 0; offset1 < 4; ++offset1) {
    for (unsigned offset2 = 0; offset2 < 4; ++offset2) {
      for (bool rc1 : {false, true}) {
        for (bool rc2 : {false, true}) {
          for (unsigned pos1 = 0; pos1 < k_max_compare_size; ++pos1) {
            dna_slice slice1 = make_slice(offset1, CHANGE_TO_SMALL, pos1, rc1);
            dna_slice slice2 = make_slice(offset2, NONE, 0, rc2);
            EXPECT_EQ(dna_compare_result::FIRST_IS_LESS,
                      slice1.compare_to(slice2))
                << "Offset1: " << offset1 << " Offset2: " << offset2
                << " Rc1: " << rc1 << " Rc2: " << rc2 << " Pos1: " << pos1
                << " Slice1:\n"
                << slice1.as_string() << " Slice2:\n"
                << slice2.as_string();
            EXPECT_EQ(dna_compare_result::SECOND_IS_LESS,
                      slice2.compare_to(slice1))
                << "Offset1: " << offset1 << " Offset2: " << offset2
                << " Rc1: " << rc1 << " Rc2: " << rc2 << " Pos1: " << pos1
                << " Slice1:\n"
                << slice1.as_string() << " Slice2:\n"
                << slice2.as_string();
            check_shared(slice1, slice2);
          }
        }
      }
    }
  }
}

TEST_F(dna_compare_test, greater_than) {
  for (unsigned offset1 = 0; offset1 < 4; ++offset1) {
    for (unsigned offset2 = 0; offset2 < 4; ++offset2) {
      for (bool rc1 : {false, true}) {
        for (bool rc2 : {false, true}) {
          for (unsigned pos1 = 0; pos1 < k_max_compare_size; ++pos1) {
            dna_slice slice1 = make_slice(offset1, CHANGE_TO_BIG, pos1, rc1);
            dna_slice slice2 = make_slice(offset2, NONE, 0, rc2);
            EXPECT_EQ(dna_compare_result::SECOND_IS_LESS,
                      slice1.compare_to(slice2))
                << "Offset1: " << offset1 << " Offset2: " << offset2
                << " Rc1: " << rc1 << " Rc2: " << rc2 << " Pos1: " << pos1
                << " Slice1:\n"
                << slice1.as_string() << " Slice2:\n"
                << slice2.as_string();
            EXPECT_EQ(dna_compare_result::FIRST_IS_LESS,
                      slice2.compare_to(slice1))
                << "Offset1: " << offset1 << " Offset2: " << offset2
                << " Rc1: " << rc1 << " Rc2: " << rc2 << " Pos1: " << pos1
                << " Slice1:\n"
                << slice1.as_string() << " Slice2:\n"
                << slice2.as_string();
            check_shared(slice1, slice2);
          }
        }
      }
    }
  }
}

TEST_F(dna_compare_test, first_is_prefix) {
  for (unsigned offset1 = 0; offset1 < 4; ++offset1) {
    for (unsigned offset2 = 0; offset2 < 4; ++offset2) {
      for (bool rc1 : {false, true}) {
        for (bool rc2 : {false, true}) {
          for (unsigned pos1 = 0; pos1 < k_max_compare_size; ++pos1) {
            dna_slice slice1 = make_slice(offset1, TRUNCATE, pos1, rc1);
            dna_slice slice2 = make_slice(offset2, NONE, 0, rc2);
            EXPECT_EQ(dna_compare_result::FIRST_IS_PREFIX,
                      slice1.compare_to(slice2))
                << "Offset1: " << offset1 << " Offset2: " << offset2
                << " Rc1: " << rc1 << " Rc2: " << rc2 << " Pos1: " << pos1
                << " Slice1:\n"
                << slice1.as_string() << " Slice2:\n"
                << slice2.as_string();
            check_shared(slice1, slice2);
          }
        }
      }
    }
  }
}

TEST_F(dna_compare_test, second_is_prefix) {
  for (unsigned offset1 = 0; offset1 < 4; ++offset1) {
    for (unsigned offset2 = 0; offset2 < 4; ++offset2) {
      for (bool rc1 : {false, true}) {
        for (bool rc2 : {false, true}) {
          for (unsigned pos2 = 0; pos2 < k_max_compare_size; ++pos2) {
            dna_slice slice1 = make_slice(offset1, NONE, 0, rc1);
            dna_slice slice2 = make_slice(offset2, TRUNCATE, pos2, rc2);
            EXPECT_EQ(dna_compare_result::SECOND_IS_PREFIX,
                      slice1.compare_to(slice2))
                << "Offset1: " << offset1 << " Offset2: " << offset2
                << " Rc1: " << rc1 << " Rc2: " << rc2 << " Pos2: " << pos2
                << " Slice1:\n"
                << slice1.as_string() << " Slice2:\n"
                << slice2.as_string();
            check_shared(slice1, slice2);
          }
        }
      }
    }
  }
}

class dna_copy_test : public testing::Test {
 protected:
  static constexpr size_t k_max_copy_size = 28 * 5;

  void for_all_slices(const std::function<void(dna_slice)>& f) {
    dna_sequence orig_seq = make_random_sequence(k_max_copy_size);
    for (size_t offset = 0; offset < 4; ++offset) {
      SCOPED_TRACE("offset: " + std::to_string(offset));
      for (size_t len = 0; len + offset < k_max_copy_size; ++len) {
        SCOPED_TRACE("length: " + std::to_string(len));
        for (bool rc : {false, true}) {
          SCOPED_TRACE(rc ? "revcomp" : "forward");
          dna_slice slice = dna_slice(orig_seq).subseq(offset, len);
          if (rc) {
            slice = slice.rev_comp();
          }
          f(slice);
        }
      }
    }
  }
};

constexpr size_t dna_copy_test::k_max_copy_size;

TEST_F(dna_copy_test, copy) {
  for_all_slices([](dna_slice slice) {
    dna_sequence result(slice);
    EXPECT_EQ(dna_compare_result::EQUAL, slice.compare_to(result)) << "slice: " << slice
                                                                   << " result: " << result;
  });
}

TEST_F(dna_copy_test, assign) {
  for_all_slices([](dna_slice slice) {
    dna_sequence result;
    result = slice;
    EXPECT_EQ(dna_compare_result::EQUAL, slice.compare_to(result))
        << "slice: " << slice << " result: " << result;
  });
}

TEST_F(dna_copy_test, append) {
  for_all_slices([](dna_slice slice) {
    for (unsigned i = 0; i < slice.size(); ++i) {
      dna_slice slice1 = slice.subseq(0, i);
      dna_slice slice2 = slice.subseq(i, slice.size() - i);
      dna_sequence result;
      result += slice1;
      EXPECT_EQ(slice1, result);
      EXPECT_EQ(dna_compare_result::EQUAL, slice1.compare_to(result)) << "slice1: " << slice1
                                                                      << " result: " << result;
      result += slice2;
      EXPECT_EQ(dna_compare_result::EQUAL, slice.compare_to(result)) << "slice: " << slice
                                                                     << " result: " << result;
    }
  });
}
