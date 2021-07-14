#pragma once

#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/web/chunked_encoding.h"
#include "modules/web/httpserver.h"
#include "modules/web/jsontypes.h"
#include "modules/web/rest_exceptions.h"

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

// Very generic rest handler
class rest_handler {
 public:
  rest_handler(http_request& request);
  virtual ~rest_handler() = default;

  virtual bool auth() { return true; }  // default: allowed

  virtual void get() { throw method_not_allowed("GET"); }
  virtual void put() { throw method_not_allowed("PUT"); }
  virtual void post() { throw method_not_allowed("POST"); }
  virtual void del() { throw method_not_allowed("DELETE"); }
  virtual void patch() { throw method_not_allowed("PATCH"); }
  virtual void options() { throw method_not_allowed("OPTIONS"); }

  // To be used by the implementors of get/put/etc
  std::string get_content_type() { return m_request.get_header("Content-Type"); }
  std::string get_match_result(int i) { return m_request.uri_match()[i]; }
  std::string get_cookie(const std::string& name) const;  // Returns empty string on no-such-cookie
  readable& get_input() { return m_request.conn(); }
  void set_status_code(int status,
                       const std::string& message = "OK");  // Not required, defaults to 200
  void set_cookie(const std::string& name, const std::string& value,
                  int max_age);  // Lets you set a cookie
  void set_content_length(
      size_t length);  // Not required, but useful to allow persistant connection

  // Allows output, throws after this are MISHANDLED, beware
  writable& set_output(const std::string& content_type, const std::string& filename = "",
                       bool use_chunked_transfer_encoding = false);

  void write_output(const std::string& content, const std::string& content_type);

 protected:
  http_request& m_request;

 private:
  chunked_encoding_writable m_chunker;
};

// Even easier rest_hander, presumes a fixed content type, defaults to application/json
class easy_rest_handler : public rest_handler {
 public:
  easy_rest_handler(http_request& request) : rest_handler(request) {}

  virtual size_t max_size() { return 64 * 1024 * 1024; }
  virtual std::string content_type() { return JSONTYPE_FULL; }

  void get() override;
  void put() override;
  void post() override;
  void del() override;
  void patch() override;
  void options() override;

  virtual std::string easy_get() { throw method_not_allowed("GET"); }
  virtual bool easy_put(const std::string& input) { throw method_not_allowed("PUT"); }
  virtual std::string easy_post(const std::string& input) { throw method_not_allowed("POST"); }
  virtual bool easy_del() { throw method_not_allowed("DELETE"); }
  virtual std::string easy_patch(const std::string& input) { throw method_not_allowed("PATCH"); }
  virtual std::string easy_options(const std::string& input) { throw method_not_allowed("OPTIONS"); }
};

std::string read_entity(http_request& request, size_t max_size);

// Called from main after all the static 'registers' are done
void run_restful_server(const bind_list_t& bind_list, const std::string& pem_path,
                        const std::string& ssl_certificates_chain_path,
                        const std::string& fork_mode, bool block = true);

typedef std::function<rest_handler*(http_request&)> create_handler_t;
void register_handler(const std::string& path, create_handler_t);

#define REST_REGISTER(classname)                                                          \
  class register_##classname {                                                            \
   public:                                                                                \
    register_##classname() { register_handler(classname::pattern(), create); }            \
    static rest_handler* create(http_request& request) { return new classname(request); } \
  } g_register_##classname

template <class T>
rest_handler* fwrap(http_request& request) {
  return new T(request);
}
