#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

std::vector<assembly_ptr> make_ref_assemblies(const scaffold& s, aoffset_t start_offset,
                                              aoffset_t end_offset, aoffset_t max_chunk_size);

}  // namespace variants
