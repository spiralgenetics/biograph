
#pragma once
#include "modules/io/version.h"

struct spec_header
{
	struct scaffold_t {
		TRANSFER_OBJECT {
			VERSION(0);
			FIELD(name);
			FIELD(md5);
			FIELD(size);
		}

		std::string name;
		std::string md5;
		size_t size;
	};

	product_version version;
	std::vector<scaffold_t> scaffolds;

	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(version);
		FIELD(scaffolds);
	}
};

struct spec_block_ref
{
	spec_block_ref() : id(0), offset(0) {}
	uint64_t id = std::numeric_limits<uint64_t>::max();
	uint64_t offset = std::numeric_limits<uint64_t>::max();
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(id);
		FIELD(offset);
	}
};

struct spec_toc
{
	spec_block_ref spec_header;
	spec_block_ref bam_header;
	spec_block_ref data_start;
	spec_block_ref index;
	spec_block_ref ref_index;
	spec_block_ref ref_data;
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(spec_header);
		FIELD(bam_header);
		FIELD(data_start);
		FIELD(index);
		FIELD(ref_index);
		FIELD(ref_data);
	}
};

struct back_ptr
{
	uint64_t id;
	uint64_t offset;
};

