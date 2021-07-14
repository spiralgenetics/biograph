#pragma once

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <cstdint>
#include <iostream>
#include <map>

// Autostats provides a way to track multiple types of metrics,
class autostats_base {
 public:
  virtual std::map<std::string, size_t> value_map() const = 0;
  void write_to_stream(std::ostream& os) const;
};

class autostats_max_value {
 public:
  void add(size_t val) { m_value = std::max(val, m_value); }
  size_t value() const { return m_value; }

 private:
  size_t m_value = 0;
};

#define DECLARE_AUTOSTATS(NAME, FIELDS)                                                       \
  NAME& operator+=(const NAME& autostats_rhs) {                                               \
    BOOST_PP_SEQ_FOR_EACH(AUTOSTATS_APPLY, ADD_FIELD, FIELDS)                                 \
    return *this;                                                                             \
  }                                                                                           \
  std::map<std::string, size_t> value_map() const override {                                  \
    std::map<std::string, size_t> autostats_result;                                           \
    BOOST_PP_SEQ_FOR_EACH(AUTOSTATS_APPLY, POPULATE_VALUE_MAP, FIELDS);                       \
    return autostats_result;                                                                  \
  }                                                                                           \
  friend std::ostream& operator<<(std::ostream& os, const NAME& st) __attribute__((unused)) { \
    st.write_to_stream(os);                                                                   \
    return os;                                                                                \
  }                                                                                           \
  BOOST_PP_SEQ_FOR_EACH(AUTOSTATS_APPLY, DECLARE_FIELD, FIELDS);

#define AUTOSTATS_APPLY(R, OPERATION, ELEM) \
  AUTOSTATS_APPLY2(OPERATION, BOOST_PP_TUPLE_ELEM(2, 0, ELEM), BOOST_PP_TUPLE_ELEM(2, 1, ELEM))
#define AUTOSTATS_APPLY2(OPERATION, COUNTER_TYPE, FIELD_NAME)                      \
  BOOST_PP_CAT(AUTOSTATS_, BOOST_PP_CAT(BOOST_PP_CAT(OPERATION, _), COUNTER_TYPE)) \
  (FIELD_NAME)

#define AUTOSTATS_DECLARE_FIELD_COUNTER(FIELD_NAME) size_t FIELD_NAME = 0;
#define AUTOSTATS_DECLARE_FIELD_MAX(FIELD_NAME) autostats_max_value FIELD_NAME;

#define AUTOSTATS_ADD_FIELD_COUNTER(FIELD_NAME) this->FIELD_NAME += autostats_rhs.FIELD_NAME;
#define AUTOSTATS_ADD_FIELD_MAX(FIELD_NAME) this->FIELD_NAME.add(autostats_rhs.FIELD_NAME.value());

#define AUTOSTATS_POPULATE_VALUE_MAP_MAX(FIELD_NAME) \
  autostats_result.insert(std::make_pair(#FIELD_NAME, this->FIELD_NAME.value()));
#define AUTOSTATS_POPULATE_VALUE_MAP_COUNTER(FIELD_NAME) \
  autostats_result.insert(std::make_pair(#FIELD_NAME, this->FIELD_NAME));
