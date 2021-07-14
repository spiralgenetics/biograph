
#include <gtest/gtest.h>
#include <map>
#include "base/base.h"
#include "modules/io/range_coder.h"
#include "modules/io/mem_io.h"
#include "modules/io/log.h"
#include "modules/io/dynamic_codecs.h"
#include "modules/io/dict_codec.h"
#include "modules/io/file_io.h"

TEST(range_coder, round_trip) {
	prefix_sum_dist dist(3);
	std::vector<double> probs = { 3. / 6., 2. / 6., 1. / 6.};
	dist.inner().add(0, 3);
	dist.inner().add(1, 2);
	dist.inner().add(2, 1);
	size_t count = 100000;
	//size_t count = 10;

	mem_io coded("", track_alloc("range_coder_test"));
	std::vector<uint8_t> actual;
	range_encoder re(coded);
	std::vector<size_t> hist(3);
	for(size_t i = 0; i < count; i++) {
		uint8_t rn = random()%6;
		rn = (rn < 3 ? 0 : (rn < 5 ? 1 : 2));
		hist[rn]++;
		re.encode(dist, rn);
		actual.push_back(rn);
	}
	re.end();
	double entropy_per_entry = 0.0;
	for(size_t i = 0; i < probs.size(); i++) {
		entropy_per_entry += -probs[i] * log2(probs[i]);
	}
	double entropy = entropy_per_entry * count;
	double size_bytes = entropy / 8.0;
	SPLOG("Expected size: %f, Coded size: %zu", size_bytes, coded.size());
	CHECK(fabs(size_bytes - double(coded.size())) < 100);
	range_decoder rd(coded);
	for(size_t i = 0; i < count; i++) {
		uint8_t rn = (uint8_t) rd.decode(dist);
		//SPLOG("%d, %d", rn, actual[i]);
		CHECK_EQ(rn, actual[i]);
	}
}

class prob_dist 
{
public:
        prob_dist(uint32_t small, uint32_t tot) : m_small(small), m_tot(tot) {}
        symbol_t count() const { return 2; }
        void symbol_range(symbol_t s, range_t r, range_t& start, range_t& end) const
        {
		range_t m = uint64_t(r) * uint64_t(m_small) / uint64_t(m_tot);
		if (s == 0) {
			start = 0;
			end = m;
		} else {
			start = m;
			end = r;
		}
        }
        void symbol_find(range_t x, range_t r, symbol_t& s, range_t& start, range_t& end) const
        {
		range_t m = uint64_t(r) * uint64_t(m_small) / uint64_t(m_tot);
		if (x < m) {
			s = 0;
		} else {
			s = 1;
		}
                symbol_range(s, r, start, end);
        }
private:
	range_t m_small;
	range_t m_tot;
};


TEST(range_coder, perf_bit) {
	size_t m100 = 10*1024*1024;
	prob_dist pd(20, 100);
	mem_io coded("", track_alloc("range_coder_test"));
	coded.reserve(m100/8);
	time_t start = time(0);
	SPLOG("Encoding 10M bits");
	range_encoder re(coded);
	for(size_t i = 0; i < m100; i++)
	{
		re.encode(pd, symbol_t(rand() % 100 < 20));
	}
	re.end();
	SPLOG("Time = %ld", (time(0) - start));
	start = time(0);
	SPLOG("Decoding 10M bits");
	range_decoder rd(coded);
	for(size_t i = 0; i < m100; i++)
	{
		symbol_t d = rd.decode(pd);
		ASSERT_LT(d, 2);
	}
	SPLOG("Time = %ld", (time(0) - start));
}

TEST(range_coder, perf) {
	size_t m100 = 100*1024*1024;
	uniform_dist u(256);
	mem_io coded("", track_alloc("range_coder_test"));
	coded.reserve(m100 + m100/20);
	time_t start = time(0);
	SPLOG("Encoding 100M bytes");
	range_encoder re(coded);
	for(size_t i = 0; i < m100; i++)
	{
		re.encode(u, i % 256);
	}
	re.end();
	SPLOG("Time = %ld", (time(0) - start));
	start = time(0);
	SPLOG("Decoding 100M bytes");
	range_decoder rd(coded);
	for(size_t i = 0; i < m100; i++)
	{
		symbol_t d = rd.decode(u);
		ASSERT_EQ(d, i % 256);
	}
	SPLOG("Time = %ld", (time(0) - start));
}

TEST(range_coder, dyn) {
	size_t count = 100000;
	mem_io coded("", track_alloc("range_coder_test"));
	dyn_prob_codec the_dp_codec(20);
	range_encoder the_range_encoder(coded);
	std::vector<symbol_t> remember;
	std::map<symbol_t, unsigned int> hist;
	SPLOG("Encoding dyn");
	time_t start = time(0);
	for(size_t i = 0; i < count; i++) {
		symbol_t the_symbol = random() % 19 >> (random() % 3);
		remember.push_back(the_symbol);
		hist[the_symbol]++;
		the_dp_codec.encode(the_range_encoder, the_symbol);
	}
	SPLOG("Encode time = %ld", (time(0) - start));
	the_range_encoder.end();
	SPLOG("Decoding dyn");
	dyn_prob_codec second_dp_codec(20);
	range_decoder the_range_decoder(coded);
	start = time(0);
	for(size_t i = 0; i < count; i++) {
		symbol_t the_symbol;
		second_dp_codec.decode(the_range_decoder, the_symbol);
		ASSERT_EQ(the_symbol, remember[i]);
	}
	SPLOG("Decode time = %ld", (time(0) - start));
	double entropy = 0.0;
	for(const auto& kvp : hist) {
		entropy += -log2(double(kvp.second)/double(count)) * kvp.second;
	}
	SPLOG("Size = %d, theoretical = %f, diff = %f", (int) coded.size()*8, entropy, coded.size()*8 - entropy);
}

TEST(range_coder, dict) {
  mem_io coded("", track_alloc("range_coder_test"));
	dict_codec d1(14);
	range_encoder re(coded);
	std::string test_msg = slurp_file("/src/golden/e_coli_10000snp.fq");
	std::string gzipped = slurp_file("/src/golden/e_coli_10000snp.fq.gz");
	size_t size = test_msg.size();
	SPLOG("Encoding");
	for(size_t i = 0; i < size; i++) {
		d1.encode(re, uint8_t(test_msg[i]));
		if (i % 10000 == 9999) {
			d1.enc_eor(re);
		}
	}
	d1.enc_eor(re);
	re.end();
	SPLOG("Original size = %d, gzip size= %d, compressed size = %d",
		(int) test_msg.size(), (int) gzipped.size(), (int) coded.size());
	dict_codec d2(14);
	range_decoder rd(coded);
	for(size_t i = 0; i < size; i++) {
		uint8_t c;
		d2.decode(rd, c);
		if (i % 10000 == 9999) {
			d2.dec_eor(rd);
		}
		if (test_msg[i] != char(c)) {
			SPLOG("Mismatch, i = %d", (int) i);
		}
		ASSERT_EQ(test_msg[i], char(c));
	}
	// ASSERT_LT(coded.size(), gzipped.size());
}

