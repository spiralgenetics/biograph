

#include "modules/bio_base/dna_sequence.h"
#include "modules/mapred/sorter.h"
#include "modules/io/log.h"

class prefix_sorter : public sorter
{
public:
        prefix_sorter(const std::string& params) {}
        int compare(const std::string& key1, const std::string& key2) const override;
	size_t partition(const std::string& key, size_t num_partitions) const override;
};

REGISTER_1(sorter, prefix, const std::string&);

int prefix_sorter::compare(const std::string& key1, const std::string& key2) const
{
	dna_sequence k1, k2;
	if (key1 != "") msgpack_deserialize(k1, key1);
	if (key2 != "") msgpack_deserialize(k2, key2);
	int r = 0;
	dna_slice s1(k1);
	dna_slice s2(k2);
	size_t msize = std::min(s1.size(), s2.size());
	bool is_sub = (s1.subseq(0, msize) == s2.subseq(0, msize));
	bool is_lt = (s1 < s2);
	if (is_lt) {
		r = is_sub ? -1 : -2;
	} else {
		r = is_sub ? (s1.size() == s2.size() ? 0 : 1) : 2;
	}
	//SPLOG("prefix_sort: k1 = '%s', k2 = '%s', r = %d", k1.as_string().c_str(), k2.as_string().c_str(), r);
	return r;
}

size_t prefix_sorter::partition(const std::string& key, size_t num_partitions) const
{
	if (num_partitions == 1) return 0;
	throw io_exception("It's invalid to partition using the prefix_sorter");
}

