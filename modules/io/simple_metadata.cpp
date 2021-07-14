#include "modules/io/simple_metadata.h"

class discard_simple_metadata_impl : public simple_metadata {
 public:
  void set_simple_json(const std::string& key, const js::mValue& value) override {}
};

simple_metadata& discard_simple_metadata() {
  static discard_simple_metadata_impl impl;
  return impl;
}

