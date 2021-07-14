#include "modules/io/encoding.h"
#include "modules/io/make_unique.h"
#include "modules/io/pass_thru.h"
#include "modules/io/zip.h"
#include "modules/io/bzip.h"

#include <map> 

const char* codec::null  = "null";
const char* codec::gzip  = "gzip";
const char* codec::bzip2 = "bzip2";

typedef std::function<std::unique_ptr<writable>(writable&)> encoder_f;
typedef std::function<std::unique_ptr<readable>(readable&)> decoder_f;

static 
const std::map<std::string, encoder_f> g_encoders {
	{ codec::null, [](writable& in) { return make_unique<pass_thru_writable>(in); } },
	{ codec::gzip, [](writable& in) { return make_unique<zip_writer>(in, no_update, 1); } },
	{ "gzip1", [](writable& in) { return make_unique<zip_writer>(in, no_update, 1); } },
	{ "gzip9", [](writable& in) { return make_unique<zip_writer>(in, no_update, 9); } },
};

static 
const std::map<std::string, decoder_f> g_decoders {
	{ "",           [](readable& in) { return make_unique<pass_thru_readable>(in); } },
	{ codec::null,  [](readable& in) { return make_unique<pass_thru_readable>(in); } },
	{ codec::gzip,  [](readable& in) { return make_unique<zip_reader>(in); } },
	{ "gzip1",  [](readable& in) { return make_unique<zip_reader>(in); } },
	{ "gzip9",  [](readable& in) { return make_unique<zip_reader>(in); } },
	{ codec::bzip2, [](readable& in) { return make_unique<bzip_reader>(in); } },
};

template <typename T, typename C> 
T find(const C& table, const std::string& encoding)
{
	const auto it = table.find(encoding);
	if (it == table.cend()) {
		throw unknown_codec(encoding);
	}
	return it->second;
}

std::unique_ptr<writable> make_encoder(const std::string& encoding, writable& sink)
{
	return find<encoder_f>(g_encoders, encoding)(sink);
}

std::unique_ptr<readable> make_decoder(const std::string& encoding, readable& source)
{
	return find<decoder_f>(g_decoders, encoding)(source);
}
