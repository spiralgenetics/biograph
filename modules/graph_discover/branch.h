#pragma once

#include <boost/optional.hpp>
#include "absl/container/btree_map.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/graph_discover/discover.h"
#include "modules/variants/assemble.h"
#include "modules/variants/ref_map.h"
#include "absl/container/flat_hash_set.h"
#include "modules/variants/scaffold.h"

namespace variants {

class branch_discover : public graph_discover {
 public:
  branch_discover(const assemble_options& options, const std::string& tag, pipeline_step_t output);
  ~branch_discover() override;

  void on_trace(const active_assembly*) override;

 private:
  std::string m_tag;
};

}  // namespace variants
