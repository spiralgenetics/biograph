#include "modules/io/config.h"
#include "modules/io/utils.h"

#include <boost/format.hpp>
#include <boost/signals2.hpp>
#include <fstream>
#include <iostream>

#define COMMENT_CHAR '#'
#define END_OF_COMMENT_CHAR '\n'

using boost::format;
using boost::str;

unknown_key_exception::unknown_key_exception(const std::string& key)
    : io_exception(str(format("Cannot find key '%s' in config") % key)) {}

Config& Config::instance() {
  static Config c;
  return c;
}

void setup_synthetic(Config& config);

void Config::load(const std::string& configfile) {

  instance().m_config_file = configfile;
  // slurp the file into a 'config' string
  std::ifstream ifs(configfile);
  if (!ifs.is_open()) {
    throw io_exception(str(format("'%s' could not be opened.") % configfile));
  }

  std::string config((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

  // remove all comments
  size_t comment_pos = config.find_first_of(COMMENT_CHAR);
  while (comment_pos != std::string::npos) {
    size_t end_of_comments_pos = config.find_first_of(END_OF_COMMENT_CHAR, comment_pos);
    config.erase(comment_pos, end_of_comments_pos - comment_pos + 1);
    comment_pos = config.find_first_of(COMMENT_CHAR);
  }

  // feed 'config' into the json_spiral reader.
  js::mValue v;
  if (!js::read(config, v)) {
    throw io_exception(str(format("'%s' is an invalid configuration file.") % configfile));
  }

  // merge the file into the existing config (where defaults may already exist)
  for(auto k : v.get_obj()) {
    instance().m_config[k.first] = k.second;
  }

  // set up dependent defaults
  setup_synthetic(instance());
}

template <class T>
static std::pair<T, bool> get(Config& config, const std::string& name) {
  const auto it = config.m_config.find(name);
  if (it == config.m_config.end()) {
    return std::make_pair(T(), false);
  }

  T value;
  json_unwrap(value, it->second);
  return std::make_pair(value, true);
}

template <class T>
static T get(Config& config, const std::string& name, const T& default_value) {
  const auto it = config.m_config.find(name);
  if (it == config.m_config.end()) {
    return default_value;
  }

  T value;
  json_unwrap(value, it->second);
  return value;
}

const char* getenv_raw(const char* name) {
  const char* var = std::getenv(name);
  if (var == nullptr) {
    throw io_exception(str(format("Missing environment variable: %s") % name));
  }
  return var;
}

std::string getenv_str(const char* name) { return std::string(getenv_raw(name)); }

std::string getenv_str(const char* name, const std::string& default_value) {
  const char* var = std::getenv(name);
  if (var == nullptr) {
    return default_value;
  }
  return var;
}

int getenv_int(const char* name) { return std::stoi(getenv_raw(name)); }

int getenv_int(const char* name, const int& default_value) {
  const char* var = std::getenv(name);
  if (var == nullptr) {
    return default_value;
  }
  return std::stoi(var);
}

// Default values that depend on the presence of other config variables.
// See config.h for static defaults.
void setup_synthetic(Config& config) {

  auto install = get<std::string>(config, "install_root", "");
  if (not install.empty()) {
    Config::set("pem_file", str(format("%s/etc/keys/ssl.pem") % install));
    Config::set("ssl_certificates_chain", str(format("%s/etc/keys/ssl.crt") % install));
  }

  auto storage = get<std::string>(config, "storage_root", "");
  if (storage.empty()) {
    storage = getenv_str("STORAGE_ROOT");
    Config::set("storage_root", storage);
  }
  auto gc_root = get<std::string>(config, "gc_root", "");
  if (gc_root.empty()) {
    Config::set("gc_root", str(format("%s/gc") % storage));
  }

  auto resources_root = get<std::string>(config, "resources_root", "");
  if (resources_root.empty()) {
    Config::set("resources_root", str(format("%s/resources") % storage));
  }

  Config::set("path_reference_base", str(format("%s/reference") % storage));
  Config::set("path_user_base", str(format("%s/users") % storage));
  Config::set("path_bulkdata", str(format("%s/bulkdata") % storage));

  // Default local temporary directory. Set in config or override with posix TMPDIR.
  auto temp_root = get<std::string>(config, "temp_root", "");
  if (temp_root.empty()) {
    auto path_bulkdata = get<std::string>(config, "path_bulkdata");
    Config::set("temp_root", path_bulkdata.first);
  }
  std::string tmpdir = getenv_str("TMPDIR", "");
  if (not tmpdir.empty()) {
    Config::set("temp_root", tmpdir);
  }
}
