#ifndef __config_h__
#define __config_h__

// Config is a singleton class (not thread-safe as of now!) designed to hold a
// list of flat constants optionally initialized from a json file. This object
// can be overridden with Config::set().

#include <functional>
#include "modules/io/json_transfer.h"
#include "modules/io/utils.h"

namespace js = json_spirit;

class unknown_key_exception : public io_exception {
 public:
  unknown_key_exception(const std::string& key);
};

const char* getenv_raw(const char* name);
std::string getenv_str(const char* name);
std::string getenv_str(const char* name, const std::string& default_value);
int getenv_int(const char* name);
int getenv_int(const char* name, const int& default_value);

class Config {
  //____ basic singleton properties ____
 public:
  static Config& instance();

  Config() {
    // config defaults
    std::string default_cfg = R"|(
      {
        "http_server_password" : "",
        "log_http_requests"    : false,
        "log_http_traffic"     : false,
        "resource_quota_slop"  : 1073741824,
        "task_update_interval" : 10,
        "url_base"             : "/api/users/",
      }
    )|";

    js::mValue v;
    js::read(default_cfg, v);
    m_config = v.get_obj();
  }

  // Config& operator=(const Config&);
  //___ end of basic singleton properties ____

 public:
  // throws io_exception if file does not exist or contains an invalid JSON object.
  static void load(const std::string& configfile);

  // throws an exception if the input param does not exist in the configuration file.
  template <class T>
  T get(const std::string& param) {
    js::mObject::const_iterator it = m_config.find(param);
    if (it == m_config.end()) throw unknown_key_exception(param);
    T dest;
    json_unwrap(dest, it->second);
    return dest;
  }

  template <class T>
  T get(const std::string& name, const T& default_value) {
    const auto it = m_config.find(name);
    if (it == m_config.end()) {
      return default_value;
    }

    T value;
    json_unwrap(value, it->second);
    return value;
  }

  template <class T>
  static void set(const std::string& name, const T& value) {
    instance().m_config[name] = js::mValue(value);
  }

  std::string m_config_file;
  js::mObject m_config;
};

// the purpose of this class is to allow implicit type detection by the compiler when calling
// Config::get. This in turn simplifies the syntax of the CONF macro, by not requiring the
// specification of the type.
class Proxy {
  std::string m_name;

 public:
  Proxy(const std::string& param) { m_name = param; }
  template <class T>
  operator T() {
    return Config::instance().get<T>(m_name);
  }
  operator std::string() { return Config::instance().get<std::string>(m_name); }
};

#define CONF(param) Proxy(#param)
#define CONF_T(type, param) Config::instance().get<type>(#param)

// a lot of config params values are strings, so the following macros simplify a common usage for
// them.
#define CONF_S(param) CONF_T(std::string, param)
#define CONF_CS(param) CONF_T(std::string, param).c_str()

#endif
