
#include <stdio.h>
#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include <map>

#include "modules/io/tunstall.h"
#include "modules/io/log.h"

TEST(tunstall, tunstall)
{
	size_t tsize = 1 << 16;

	// Validate creation and writing
	tunstall t(.05, tsize);
	ASSERT_EQ(t.size(), tsize);
	uint8_t *buf = new uint8_t[tunstall::buf_size(tsize)];
	t.write(buf);
	tunstall t2(buf, tunstall::buf_size(tsize));
	ASSERT_EQ(t2.size(), tsize);
	for(size_t i = 0; i < t.size(); i++) {
		ASSERT_EQ(t[i], t2[i]);
	}

	// Validate encoding
	// Make zeroed out data
	std::vector<uint8_t> data(1024);
	// Or in 1/20th of the bits or so
	for(size_t i = 0; i < 1024*8/20; i++) {
		data[rand()%1024] |= (1 << (rand() % 8));
	}
	// Encode
	std::vector<uint16_t> enc_data;
	t.encode(enc_data, &data[0], 1024);
	SPLOG("Expected data size ~= %f", (-log2(.05)*.05 - log2(.95)*.95)*8192.0/8.0);
	SPLOG("Size of encoded data = %lu bytes", enc_data.size()*2);
	std::vector<uint8_t> dec_data(1024);
	t.decode(enc_data, &dec_data[0], 1024);
	for(size_t i = 0; i < 1024; i++) {
		if (data[i] != dec_data[i]) {
			SPLOG("%lu: %x vs %x", i, data[i], dec_data[i]);
		}
	}
	ASSERT_EQ(data, dec_data);
}

