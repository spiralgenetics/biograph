
#ifndef __sorter_h__
#define __sorter_h__

#include <string>

#include "modules/io/msgpack_transfer.h"
#include "modules/io/registry.h"

class sorter
{
public:
	virtual ~sorter() {}
	// Returns:
	// +-2 if different group, -2 if k1 < k2, +2 if k1 > k2
	// +-1 if same group, -1 if k1 < k2, +1 if k1 > k2
	// 0 is identical
	virtual int compare(const std::string& key1, const std::string& key2) const = 0;
	virtual std::string bump_back(const std::string& key) const { return key; }
	virtual size_t partition(const std::string& key, size_t num_partitions) const = 0;

	inline bool gt(const std::string& key1, const std::string& key2) { return this->compare(key1, key2) >= 1; }
	inline bool lt(const std::string& key1, const std::string& key2) { return this->compare(key1, key2) <= -1; }
};

template<class T>
class simple_sorter : public sorter
{
public:
	simple_sorter(const std::string& params) {}
	int compare(const std::string& key1, const std::string& key2) const override
	{
		T a;
		T b;
		msgpack_deserialize(a, key1);
		msgpack_deserialize(b, key2);	
		if (a < b) return -2;
		if (b < a) return 2;
		return 0;
	}

	size_t partition(const std::string& key, size_t num_partitions) const override
	{
		return 0; // Correct but terrible is using partitioning
	}
};

DECLARE_REGISTRY_1(sorter, std::string const&);

#endif
