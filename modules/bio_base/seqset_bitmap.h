
#pragma once

#include "modules/io/tunstall.h"
#include "modules/io/mmap_buffer.h"

class seqset_bitmap_base
{
public:
	virtual ~seqset_bitmap_base() = default;
	virtual bool get_bit(uint64_t loc) const = 0;
};

class seqset_bitmap_true : public seqset_bitmap_base
{
public:
	seqset_bitmap_true() {}
	bool get_bit(uint64_t loc) const override { return true; }
};

