#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/make_unique.h"

// DECLARE_REGISTRY_n declares a new registry that returns the given
// type, created with the specified number of parameters.
//
// For instance, REGISTRY_1(foo, bar) declares a new registry named
// "foo_registry" that creates objects of type "foo" from subclasses
// that take a single argument, "bar".
//
// You can create an instance of subclass "baz" by calling
// foo_registry::get("baz",my_bar).  Calling get_safe instead of get
// will result in an exception being thrown if there's no "baz"
// subclass registered.
//
// When declaring a new registry, you will also need to add an
// definition for the registry using DEFINE_REGISTRY.
#define DECLARE_REGISTRY_0(name)        \
  extern template class registry<name>; \
  using name##_registry = registry<name>
#define DECLARE_REGISTRY_1(name, param1)        \
  extern template class registry<name, param1>; \
  using name##_registry = registry<name, param1>
#define DECLARE_REGISTRY_2(name, param1, param2)        \
  extern template class registry<name, param1, param2>; \
  using name##_registry = registry<name, param1, param2>
#define DECLARE_REGISTRY_3(name, param1, param2, param3)        \
  extern template class registry<name, param1, param2, param3>; \
  using name##_registry = registry<name, param1, param2, param3>

#define DEFINE_REGISTRY_0(name) template class registry<name>;
#define DEFINE_REGISTRY_1(name, param1) template class registry<name, param1>;
#define DEFINE_REGISTRY_2(name, param1, param2) template class registry<name, param1, param2>;
#define DEFINE_REGISTRY_3(name, param1, param2, param3) \
  template class registry<name, param1, param2, param3>;

// Registers a subclass creator.  For instance, REGISTER_1(foo, baz,
// bar) registers a subclass of "foo" named "baz" that can be created
// with an argument named "bar".
#define REGISTER_0(basename, name) \
  static register_helper<basename, name##_##basename> __reg_##name##basename(#name)
#define REGISTER_1(basename, name, param1) \
  static register_helper<basename, name##_##basename, param1> __reg_##name##basename(#name)
#define REGISTER_2(basename, name, param1, param2) \
  static register_helper<basename, name##_##basename, param1, param2> __reg_##name##basename(#name)
#define REGISTER_3(basename, name, param1, param2, param3)                    \
  static register_helper<basename, name##_##basename, param1, param2, param3> \
      __reg_##name##basename(#name)

// Implementation follows.
class register_helper_base {
 public:
  constexpr register_helper_base(const char* name) : m_name(name) {}
  register_helper_base* next() const { return m_next; }
  virtual ~register_helper_base();

  void set_next(register_helper_base* next) {
    CHECK(!m_next);
    m_next = next;
  }
  const char* name() const { return m_name; }

 private:
  const char* m_name = 0;
  register_helper_base* m_next = nullptr;
};

template <typename rtype, typename... param>
class registry {
 public:
  static std::unique_ptr<rtype> get(const std::string& name, param... p);
  static std::unique_ptr<rtype> get_safe(const std::string& name, param... p);

  static void add(register_helper_base* helper);

 private:
  static register_helper_base* g_registered;
};

template <typename rtype, typename... param>
register_helper_base* registry<rtype, param...>::g_registered = nullptr;

template <class rtype, typename... param>
class typed_register_helper_base : public register_helper_base {
 public:
  constexpr typed_register_helper_base(const char* name) : register_helper_base(name) {}
  virtual std::unique_ptr<rtype> my_make(param... p) = 0;
};

template <class rtype, class dtype, typename... param>
class register_helper : public typed_register_helper_base<rtype, param...> {
 public:
  std::unique_ptr<rtype> my_make(const param... p) override { return make_unique<dtype>(p...); }
  register_helper(const char* name) : typed_register_helper_base<rtype, param...>(name) {
    registry<rtype, param...>::add(this);
  }
};

template <typename rtype, typename... param>
std::unique_ptr<rtype> registry<rtype, param...>::get(const std::string& name, param... p) {
  register_helper_base* helper = g_registered;
  while (helper && helper->name() != name) {
    helper = helper->next();
  }
  if (helper) {
    using typed_base = typed_register_helper_base<rtype, param...>;
    typed_base* typed_helper = dynamic_cast<typed_base*>(helper);
    CHECK(typed_helper) << typeid(*helper).name();
    return typed_helper->my_make(p...);
  } else {
    return nullptr;
  }
}

template <typename rtype, typename... param>
std::unique_ptr<rtype> registry<rtype, param...>::get_safe(const std::string& name, param... p) {
  std::unique_ptr<rtype> r = get(name, p...);
  if (r == nullptr) {
    throw io_exception("Unknown type in registry: '" + name + "'");
  }
  return r;
}

template <typename rtype, typename... param>
void registry<rtype, param...>::add(register_helper_base* helper) {
  helper->set_next(g_registered);
  g_registered = helper;
}
