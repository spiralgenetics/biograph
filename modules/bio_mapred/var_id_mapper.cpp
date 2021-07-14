#include "modules/bio_mapred/var_id_mapper.h"

REGISTER_1(mapper, var_id, const std::string&);

void var_id_mapper::typed_map(const seq_position& key, const struct_var& var)
{
	struct_var_key sk(var.var_id, 0);
	if (var.is_structural && var.ref_start > var.ref_end) return;
	output(sk, var);
}

