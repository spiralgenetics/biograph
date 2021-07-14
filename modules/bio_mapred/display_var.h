
#pragma once 

#include "modules/io/color_text_buffer.h"
#include "modules/mapred/path.h"
#include "modules/bio_base/pileup.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_base/reference.h"

void display_var(
	color_text_buffer& out, 
	const reference& ref,
	const std::vector<struct_var>& vars,
	const std::vector<read_support>& reads
);

void display_vars(const path& out_path, const reference& ref, kv_source& svs, kv_source& reads, unsigned int min_depth);

