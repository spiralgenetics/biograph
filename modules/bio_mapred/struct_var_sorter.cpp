#include "modules/io/msgpack_transfer.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_mapred/struct_var_sorter.h"

REGISTER_1(sorter, struct_var, const std::string&);

int struct_var_sorter::compare(const std::string& key1, const std::string& key2) const
{
	struct_var_key s1;
	struct_var_key s2;
	msgpack_deserialize(s1, key1);
	msgpack_deserialize(s2, key2);

	if (s1.variation_id != s2.variation_id) {
		if (s1.variation_id < s2.variation_id)
			return -2;
		else
			return 2;
	} else {
		if (s1.read_id < s2.read_id)
			return -1;
		else if (s1.read_id > s2.read_id)
			return 1;
		else
			return 0;
	}
}

size_t struct_var_sorter::partition(const std::string& key, size_t num_partitions) const
{
	if (num_partitions == 1) return 0;
	std::string error_string = "struct_var_sorter encountered a manifest with partition count ";
	error_string += std::to_string(num_partitions);
	throw io_exception(error_string);
}

