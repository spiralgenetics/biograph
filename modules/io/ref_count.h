#pragma once

#include "base/base.h"

// Move-only reference counted item.
//
// Users must explicitly do any incrementing or decrementing.  This wrapper enforces that.
template <typename T, bool k_atomic, bool k_allow_implicit_copy, bool k_allow_implicit_delete>
class explicit_shared_ptr {
 public:
  using element_type = T;

  explicit_shared_ptr() = default;
  explicit_shared_ptr(std::nullptr_t) {}
  explicit_shared_ptr(std::unique_ptr<T> rhs) : m_count(new size_t{1}), m_elem(rhs.release()) {}
  explicit_shared_ptr(const explicit_shared_ptr& rhs) {
    static_assert(k_allow_implicit_copy, "Implicit copy disallowed");
    (*this) = rhs.clone();
  }
  explicit explicit_shared_ptr(T*) {
    // This constructor is required for explicit_shared_ptr to comply
    // with pybind11's "pointer holder" API.  However, we never use it
    // and always explicitly create the shared pointers manually.
    static_assert(k_allow_implicit_copy, "Implicit copy disallowed");
    LOG(FATAL)
        << "explicit_shared_ptr should not be constructed from an ambiguously owned pointer.";
  }
  explicit_shared_ptr(explicit_shared_ptr&& rhs) {
    m_count = rhs.m_count;
    m_elem = rhs.m_elem;
    rhs.m_count = nullptr;
    rhs.m_elem = nullptr;
  }
  explicit_shared_ptr& operator=(const explicit_shared_ptr& rhs) {
    static_assert(k_allow_implicit_copy, "Implicit copy disallowed");
    (*this) = rhs.clone();
    return (*this);
  }
  explicit_shared_ptr& operator=(explicit_shared_ptr&& rhs) {
    if (!k_allow_implicit_delete) {
      CHECK(!m_elem) << "Must explicitly discard old pointer before overwriting";
    }
    if (m_elem) {
      release_and_discard();
    }
    m_count = rhs.m_count;
    m_elem = rhs.m_elem;
    rhs.m_count = nullptr;
    rhs.m_elem = nullptr;
    return *this;
  }
  ~explicit_shared_ptr() {
    if (!k_allow_implicit_delete) {
      CHECK(!m_elem) << "Element should be released when done.";
    }
    if (m_elem) {
      release_and_discard();
    }
  }

  template <typename... Args>
  static explicit_shared_ptr make_shared(Args&&... args) {
    explicit_shared_ptr result(std::unique_ptr<T>(new T(std::forward<Args>(args)...)));
    return result;
  }

  void release_and_discard() { release().reset(); }

  std::unique_ptr<T> release() WARN_UNUSED {
    std::unique_ptr<T> result;
    CHECK(m_count);
    CHECK(m_elem);

    CHECK_GT(*m_count, 0);
    if (decref() == 0) {
      delete m_count;
      result.reset(m_elem);
    }
    m_count = nullptr;
    m_elem = nullptr;

    return result;
  }

  explicit_shared_ptr clone() const {
    CHECK(m_count);
    CHECK(m_elem);

    incref();

    explicit_shared_ptr cloned;
    cloned.m_count = m_count;
    cloned.m_elem = m_elem;

    return cloned;
  }

  T* operator->() const { return m_elem; }
  T& operator*() const { return *m_elem; }
  T* get() const { return m_elem; }
  explicit operator bool() const { return m_elem; }

  size_t use_count() const {
    CHECK(m_count);
    return *m_count;
  }

  bool operator==(const explicit_shared_ptr& rhs) const { return get() == rhs.get(); }
  bool operator!=(const explicit_shared_ptr& rhs) const { return get() != rhs.get(); }
  bool operator<=(const explicit_shared_ptr& rhs) const { return get() <= rhs.get(); }
  bool operator>=(const explicit_shared_ptr& rhs) const { return get() >= rhs.get(); }
  bool operator<(const explicit_shared_ptr& rhs) const { return get() < rhs.get(); }
  bool operator>(const explicit_shared_ptr& rhs) const { return get() > rhs.get(); }
  template <typename H>
  friend H AbslHashValue(H h, const explicit_shared_ptr& c) {
    return H::combine(std::move(h), c.get());
  }

  bool operator==(const T* rhs) const { return get() == rhs; }
  bool operator!=(const T* rhs) const { return get() != rhs; }
  bool operator<=(const T* rhs) const { return get() <= rhs; }
  bool operator>=(const T* rhs) const { return get() >= rhs; }
  bool operator<(const T* rhs) const { return get() < rhs; }
  bool operator>(const T* rhs) const { return get() > rhs; }

  friend bool operator==(const T* lhs, const explicit_shared_ptr& rhs) { return lhs == rhs.get(); }
  friend bool operator!=(const T* lhs, const explicit_shared_ptr& rhs) { return lhs != rhs.get(); }
  friend bool operator<=(const T* lhs, const explicit_shared_ptr& rhs) { return lhs <= rhs.get(); }
  friend bool operator>=(const T* lhs, const explicit_shared_ptr& rhs) { return lhs >= rhs.get(); }
  friend bool operator<(const T* lhs, const explicit_shared_ptr& rhs) { return lhs < rhs.get(); }
  friend bool operator>(const T* lhs, const explicit_shared_ptr& rhs) { return lhs > rhs.get(); }

  struct less_than {
    using is_transparent = void;

    bool operator()(const explicit_shared_ptr& lhs, const explicit_shared_ptr& rhs) const {
      return lhs < rhs;
    }
    bool operator()(const T* lhs, const explicit_shared_ptr& rhs) const { return lhs < rhs; }
    bool operator()(const explicit_shared_ptr& lhs, const T* rhs) const { return lhs < rhs; }
    bool operator()(const T* lhs, const T* rhs) const { return lhs < rhs; }
  };

 private:
  void incref() const {
    CHECK(m_count);
    if (k_atomic) {
      __sync_fetch_and_add(m_count, 1);
    } else {
      ++*m_count;
    }
  }

  size_t decref() {
    CHECK(m_count);
    size_t new_val;
    CHECK_GT(*m_count, 0);
    if (k_atomic) {
      new_val = __sync_sub_and_fetch(m_count, 1);
    } else {
      new_val = --*m_count;
    }
    return new_val;
  }

  size_t* m_count = nullptr;
  T* m_elem = nullptr;
};
