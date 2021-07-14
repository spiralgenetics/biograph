#pragma once

#include "modules/io/json_transfer.h"

#include <string>

class simple_metadata {
 public:
  virtual void set_simple_json(const std::string& key, const js::mValue&) = 0;
  virtual ~simple_metadata() = default;

  template<class V>
  void set_simple(const std::string& key, const V& value) {
    set_simple_json(key, json_wrap(const_cast<V&>(value)));
  }
};

simple_metadata& discard_simple_metadata();
