#pragma once

#include "modules/io/io.h" 
#include <memory>

namespace codec
{

extern const char* null;
extern const char* gzip;
extern const char* bzip2;

} // namespace codec

// throw if 'encoding' is an unknown codec
std::unique_ptr<writable> make_encoder(const std::string& encoding, writable& sink);

// throw if 'encoding' is an unknown codec
std::unique_ptr<readable> make_decoder(const std::string& encoding, readable& source);

class unknown_codec : public io_exception
{
public:
	unknown_codec(const std::string& encoding)
		: io_exception("Unknown codec: " + encoding) 
	{}
};
