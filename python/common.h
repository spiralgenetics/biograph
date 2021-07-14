#include <pybind11/pybind11.h>
#include <boost/optional.hpp>

#include "modules/io/progress.h"

// Executes a lambda in a separate thread, while supplying progress
// updates back to python in the original thread.
//
// Caller is responsible for releasing GIL before calling.
void execute_with_py_progress(pybind11::object py_progress,
                              const std::function<void(progress_handler_t progress)> f);

// Convenience func to generate a result for __str__ from operator<<(std::ostream&, T)
template <typename T>
std::string str_from_ostream(const T& val) {
  std::stringstream s;
  s << val;
  return s.str();
}

namespace pybind11 {
namespace detail {

// TODO(nils): We should use the optional caster once
// https://github.com/pybind/pybind11/issues/1919 is fixed.  Until
// then, include the fixed version here:

template <typename T>
struct pybind11_issue_1919_optional_caster {
  using value_conv = make_caster<typename T::value_type>;
  template <typename T_>
  static handle cast(T_ &&src, return_value_policy policy, handle parent) {
    if (!src) return none().inc_ref();
    using dereference_type = decltype(*std::forward<T_>(src));
    policy = return_value_policy_override<dereference_type>::policy(policy);
    return value_conv::cast(*std::forward<T_>(src), policy, parent);
  }

  bool load(handle src, bool convert) {
    if (!src) {
      return false;
    } else if (src.is_none()) {
      return true;  // default-constructed value is already empty
    }
    value_conv inner_caster;
    if (!inner_caster.load(src, convert)) return false;
    value.emplace(cast_op<typename T::value_type &&>(std::move(inner_caster)));
    return true;
  }
  PYBIND11_TYPE_CASTER(T, _("Optional[") + value_conv::name + _("]"));
};

template <typename T>
struct type_caster<boost::optional<T>> : pybind11_issue_1919_optional_caster<boost::optional<T>> {};

}  // namespace detail
}  // namespace pybind11
