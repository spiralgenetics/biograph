#include "modules/bio_base/karyotype_compat.h"

bool kt_supercontig::operator<(const kt_supercontig& rhs) const
{
	if (chr != rhs.chr)
		return chr < rhs.chr;
	return offset < rhs.offset;
}

kt_supercontig::kt_supercontig(const std::string& _chr, size_t _offset, size_t _len)
	: chr(_chr)
	, name(printstring("%s:%d", _chr.c_str(), (int) _offset))
	, offset(_offset)
	, len(_len)
{}
