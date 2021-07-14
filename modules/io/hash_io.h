#pragma once

#include "modules/io/io.h"

#include <openssl/md5.h>
#include <iomanip>
#include <sstream>

class md5_hash_writer : public writable 
{
public:
	md5_hash_writer()
	{
		reset();
	}

	void reset()
	{
		MD5_Init(&m_ctx);
	}

	void write(const char* buf, size_t len) override
	{
		MD5_Update(&m_ctx, buf, len);
	}

	void finish()
	{		
		MD5_Final(m_digest.data(), &m_ctx);
	}

	std::array<uint8_t, MD5_DIGEST_LENGTH> digest() const
	{
		return m_digest;
	}

	std::string hex() const
	{
		std::stringstream ss;
		ss << std::hex << std::setfill('0');
		for (const auto& octet : m_digest) {
			ss << std::setw(2) << static_cast<int>(octet);
		}
		return ss.str();

	}

	std::string etag() const
	{
		std::stringstream ss;
		ss << '"' << hex() << '"';
		return ss.str();
	}

private:
	MD5_CTX m_ctx;
	std::array<uint8_t, MD5_DIGEST_LENGTH> m_digest;
};
