#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "modules/variants/assemble.h"

// tempalted version of assembly_ptr as to not confuse the PYBIND11_DECLARE_HOLDER_TYPE macro.
template <typename T>
using assembly_ptr_tmpl = explicit_shared_ptr<T, true, true, true>;

PYBIND11_DECLARE_HOLDER_TYPE(T, ::assembly_ptr_tmpl<T>);

// Wrapper to hold dict handle for assemblies that makes sure to
// acquire the GIL whenever it increments or decrements the refcount.
class assembly_dict_holder {
 public:
  assembly_dict_holder() = delete;
  assembly_dict_holder(pybind11::dict d) { m_py_obj = d.release().ptr(); }
  assembly_dict_holder(const assembly_dict_holder& rhs) : m_py_obj(rhs.m_py_obj) {
    if (m_py_obj) {
      pybind11::gil_scoped_acquire gil;
      Py_INCREF(m_py_obj);
    }
  }
  assembly_dict_holder& operator=(const assembly_dict_holder& rhs) {
    if (m_py_obj || rhs.m_py_obj) {
      pybind11::gil_scoped_acquire gil;
      if (m_py_obj) {
        Py_DECREF(m_py_obj);
      }
      m_py_obj = rhs.m_py_obj;
      if (m_py_obj) {
        pybind11::gil_scoped_acquire gil;
        Py_INCREF(m_py_obj);
      }
    }
    return *this;
  }
  pybind11::dict get() const {
    CHECK(m_py_obj);
    return pybind11::reinterpret_borrow<pybind11::dict>(m_py_obj);
  }
  ~assembly_dict_holder() {
    if (m_py_obj) {
      pybind11::gil_scoped_acquire gil;
      Py_DECREF(m_py_obj);
    }
  }

 private:
  PyObject* m_py_obj = nullptr;
};

namespace pybind11 {
namespace detail {

// Custom assembly_ptr C++ <-> python caster to make sure __dict__ gets preserved.
template <>
class type_caster<::variants::assembly_ptr>
    : public copyable_holder_caster<::variants::assembly, ::variants::assembly_ptr> {
 public:
  using base = copyable_holder_caster<::variants::assembly, ::variants::assembly_ptr>;

  bool load(handle src, bool b) {
    bool loaded = base::load(src, b);
    if (loaded) {
      using ptr_type = ::variants::assembly_ptr&;
      ::variants::assembly_ptr& a = ptr_type(*this);
      propagate_dict(a.get(), src);
    }
    return loaded;
  }

  template <typename T>
  static handle cast(T&& src, return_value_policy policy, handle parent) {
    ::variants::assembly* aptr = src.get();
    using dereference_type = decltype(*std::forward<T>(src));
    policy = return_value_policy_override<dereference_type>::policy(policy);
    handle casted = base::cast(std::forward<T>(src), policy, parent);

    propagate_dict(aptr, casted);
    return casted;
  }

 private:
  // Make sure the python object's __dict__ and the assembly's
  // user_data are shared, so that properties are persistent after
  // assemblies go from python->c++->python.
  static void propagate_dict(::variants::assembly* aptr, handle casted) {
    if (!aptr) {
      return;
    }
    if (aptr->user_data.empty()) {
      aptr->user_data = assembly_dict_holder(casted.attr("__dict__"));
    }

    try {
      dict d = boost::any_cast<const assembly_dict_holder&>(aptr->user_data).get();
      casted.attr("__dict__") = d;
    } catch (const boost::bad_any_cast&) {
      LOG(FATAL) << "Casting an assembly ptr without an associated python dictionary";
    }
  }
};

template <>
struct type_caster<::variants::optional_aoffset> {
  using value_conv = make_caster<::variants::aoffset_t>;
  template <typename T_>
  static handle cast(T_&& src, return_value_policy policy, handle parent) {
    if (!src) return none().inc_ref();
    policy = return_value_policy_override<::variants::aoffset_t>::policy(policy);
    return value_conv::cast(::variants::aoffset_t(src), policy, parent);
  }

  bool load(handle src, bool convert) {
    if (!src) {
      return false;
    } else if (src.is_none()) {
      return true;  // default-constructed value is already empty
    }
    value_conv inner_caster;
    if (!inner_caster.load(src, convert)) return false;
    value = inner_caster;
    return true;
  }
  PYBIND11_TYPE_CASTER(::variants::optional_aoffset, _("OptionalAoffset"));
};

}  // namespace detail
}  // namespace pybind11

void bind_assembly(pybind11::module& m);
