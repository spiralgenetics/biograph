#ifndef __karyotype_compat_h__
#define __karyotype_compat_h__

#include "modules/bio_base/reference_assembly.h"
#include "modules/io/transfer_object.h"

struct kt_supercontig
{
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(chr);
		FIELD(name);
		FIELD(offset);
		FIELD(len);
	}

	std::string chr;
	std::string name;
	size_t offset;
	size_t len;

	bool operator<(const kt_supercontig& rhs) const;

	kt_supercontig() : offset(0), len(0) {} // For deserialization
	kt_supercontig(const std::string& _chr, size_t _offset, size_t _len);
};

struct kt_compat
{
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(supercontigs);
		FIELD(chromosomes);
		FIELD(chr_order);
	}

	std::set<kt_supercontig> supercontigs;
	std::set<scaffold> chromosomes;
	std::vector<std::string> chr_order;
};

#endif
