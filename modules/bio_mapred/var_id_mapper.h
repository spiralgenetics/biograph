
#ifndef __var_id_mapper_h__
#define __var_id_mapper_h__

#include "modules/mapred/mapper.h"
#include "modules/bio_base/seq_position.h"
#include "modules/bio_base/struct_var.h"

class var_id_mapper :  public typed_mapper<var_id_mapper, 
        seq_position, struct_var, 
        struct_var_key, struct_var>
{
public:
	var_id_mapper(const std::string& params) {}
	void typed_map(const seq_position& key, const struct_var& value);
};

#endif
