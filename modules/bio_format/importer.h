#pragma once

#include "modules/io/keyvalue.h"
#include "modules/io/registry.h"
#include "modules/io/simple_metadata.h"

class importer {
 public:
  virtual ~importer() = default;
  virtual void import(kv_sink& sink, simple_metadata& meta = discard_simple_metadata()) = 0;
};

DECLARE_REGISTRY_3(importer, readable&, bool, std::string const&);
