
#include "modules/io/crypto.h"
#include "modules/io/hexdump.h"
#include "modules/io/log.h"
#include <gtest/gtest.h>

const char* test_salt1 = "0123456789abcdef";
const char* test_salt2 = "0132456789abcdef";
const char* test_key1 = "Hello World";
const char* test_key2 = "Goodbye World";
const char* test_msg = "This is a test, do not pass go, do not collect $200.  This is only a test";

size_t test_msg_size = strlen(test_msg);

// These test do *NOT* test AES or GCM, which are assumed to be correct
// However they do sanity check the wrappers

// Check that round-trip works, and encrypted output is actually different then input
TEST(crypto, roundtrip)
{
	// Encrypt
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	char tag[16];
	char emesg[test_msg_size];
	SPLOG("Encrypting: %s", hexdump(std::string(test_msg, test_msg_size)).c_str());
	ctx.encrypt(23, tag, emesg, test_msg, test_msg_size);
	SPLOG("Output: %s", hexdump(std::string(emesg, test_msg_size)).c_str());
	// Validate something happened
	ASSERT_TRUE(memcmp(emesg, test_msg, test_msg_size) != 0);

	// Undo it from scratch
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	char dmesg[test_msg_size];
	bool r = ctx2.decrypt(23, tag, dmesg, emesg, test_msg_size);
	SPLOG("Decrypted: %s", hexdump(std::string(dmesg, test_msg_size)).c_str());

	// Validate it round tripped
	ASSERT_TRUE(r);
	ASSERT_TRUE(memcmp(dmesg, test_msg, test_msg_size) == 0);
}

// Check that changes for salt, key, iv, messgage or tag all break decrypt

TEST(crypto, bad_salt)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	char tag[16];
	char emesg[test_msg_size];
	ctx.encrypt(23, tag, emesg, test_msg, test_msg_size);
	crypto_ctx ctx2(test_salt2, test_key1, strlen(test_key1));
	char dmesg[test_msg_size];
	bool r = ctx2.decrypt(23, tag, dmesg, emesg, test_msg_size);
	ASSERT_FALSE(r);
}

TEST(crypto, bad_key)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	char tag[16];
	char emesg[test_msg_size];
	ctx.encrypt(23, tag, emesg, test_msg, test_msg_size);
	crypto_ctx ctx2(test_salt1, test_key2, strlen(test_key2));
	char dmesg[test_msg_size];
	bool r = ctx2.decrypt(23, tag, dmesg, emesg, test_msg_size);
	ASSERT_FALSE(r);
}

TEST(crypto, bad_id)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	char tag[16];
	char emesg[test_msg_size];
	ctx.encrypt(23, tag, emesg, test_msg, test_msg_size);
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	char dmesg[test_msg_size];
	bool r = ctx2.decrypt(24, tag, dmesg, emesg, test_msg_size);
	ASSERT_FALSE(r);
}

TEST(crypto, bad_mesg)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	char tag[16];
	char emesg[test_msg_size];
	ctx.encrypt(23, tag, emesg, test_msg, test_msg_size);
	emesg[0] ^= 1;
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	char dmesg[test_msg_size];
	bool r = ctx2.decrypt(23, tag, dmesg, emesg, test_msg_size);
	ASSERT_FALSE(r);
}

TEST(crypto, bad_tag)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	char tag[16];
	char emesg[test_msg_size];
	ctx.encrypt(23, tag, emesg, test_msg, test_msg_size);
	tag[0] ^= 1;
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	char dmesg[test_msg_size];
	bool r = ctx2.decrypt(23, tag, dmesg, emesg, test_msg_size);
	ASSERT_FALSE(r);
}

TEST(crypto, block)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	SPLOG("Encryped block:\n%s", hexdump(std::string(cblock.buffer(), cblock.size())).c_str());
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	ctx2.decrypt(nblock, 23, cblock);
	ASSERT_EQ(nblock.size(), 11);
	ASSERT_TRUE(memcmp(nblock.buffer(), "Hello world", 11) == 0);
}

TEST(crypto, block_bad_key)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt1, test_key2, strlen(test_key2));
	try {
		ctx2.decrypt(nblock, 23, cblock);
	} catch(const io_exception&) {
		return;
	}
	ASSERT_TRUE(false);
}

TEST(crypto, block_bad_salt)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt2, test_key1, strlen(test_key1));
	try {
		ctx2.decrypt(nblock, 23, cblock);
	} catch(const io_exception&) {
		return;
	}
	ASSERT_TRUE(false);
}

TEST(crypto, block_bad_iv)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	try {
		ctx2.decrypt(nblock, 24, cblock);
	} catch(const io_exception&) {
		return;
	}
	ASSERT_TRUE(false);
}

TEST(crypto, block_corrupt1)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	cblock.buffer()[10] ^= 1;
	try {
		ctx2.decrypt(nblock, 23, cblock);
	} catch(const io_exception&) {
		return;
	}
	ASSERT_TRUE(false);
}

TEST(crypto, block_corrupt2)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	cblock.buffer()[20] ^= 1;
	try {
		ctx2.decrypt(nblock, 23, cblock);
	} catch(const io_exception&) {
		return;
	}
	ASSERT_TRUE(false);
}

TEST(crypto, block_corrupt3)
{
	crypto_ctx ctx(test_salt1, test_key1, strlen(test_key1));
	mem_io block("", track_alloc("crypto_test"));
	mem_io cblock("", track_alloc("crypto_test"));
	block.print("Hello world");
	ctx.encrypt(cblock, 23, std::move(block));
	mem_io nblock("", track_alloc("crypto_test"));
	crypto_ctx ctx2(test_salt1, test_key1, strlen(test_key1));
	cblock.buffer()[30] ^= 1;
	try {
		ctx2.decrypt(nblock, 23, cblock);
	} catch(const io_exception&) {
		return;
	}
	ASSERT_TRUE(false);
}


