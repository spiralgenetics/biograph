// Generic wrappers for starting and stopping the http server,
// registering URI handlers, and handling http client requests
// Actual work is done by libpoco

#pragma once

#include <signal.h>
#include <mutex>
#include <sstream>
#include <vector>

#include <boost/regex.hpp>

#include "Poco/Base64Decoder.h"
#include "Poco/Net/HTTPRequestHandler.h"
#include "Poco/Net/HTTPRequestHandlerFactory.h"
#include "Poco/Net/HTTPResponse.h"
#include "Poco/Net/HTTPServer.h"
#include "Poco/Net/HTTPServerRequest.h"
#include "Poco/Net/HTTPServerResponse.h"
#include "Poco/Net/ServerSocket.h"
#include "Poco/Util/ServerApplication.h"

#include "modules/io/io.h"
#include "modules/io/transfer_object.h"
#include "modules/web/url_query.h"

using namespace Poco::Net;

struct bind_info {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(ip);
    FIELD(ssl, false);
    FIELD(port, (ssl ? 443 : 80));
  }
  std::string ip;  // If empty string or 0.0.0.0, IP_ANY
  bool ssl;
  int port;
};

typedef std::vector<bind_info> bind_list_t;

class connection : public read_wrapper, public write_wrapper {
 public:
  connection(HTTPServerRequest& req, HTTPServerResponse& resp) : m_req(req), m_resp(resp) {}

  void finish_headers();
  void mark_valid() { m_valid = true; }

  HTTPServerRequest& m_req;
  HTTPServerResponse& m_resp;

 private:
  friend class http_request;
  int base_read(char* buf, size_t len) override;
  int base_write(const char* buf, int len) override;
  int base_close() override { return 0; }

  bool m_valid = false;
  std::ostream* m_send = nullptr;
};

class http_request;

void write_access_header(http_request& request);
void error_response(http_request& request, unsigned code, const std::string& msg);

class handler {
 public:
  virtual void handle(http_request& request) = 0;
};

typedef std::function<bool(http_request& request)> handler_fn;

// Start and stop a libpoco http server
// Implement server request handler registry
class http_server {
  friend class http_request;

 public:
  static http_server& get();

  // Registers a handle, throw if bad regex
  void register_handler(handler& handler, const std::string& uri_regex,
                        const std::string& method_regex);

  // Star the server, bound to the ip:port(s) specified in bind_list
  void start(const bind_list_t& bind_list = bind_list_t(), const std::string& pem_path = "",
             const std::string& ssl_certificates_chain_path = "");

  void stop();

 private:
  void handle(http_request& req);
  friend class RequestHandler;

  struct handler_reg {
    handler_reg(handler& handler, const std::string& uri_regex, const std::string& method_regex);

    bool match(http_request& req);

    handler& m_handler;
    boost::regex uri;
    boost::regex method;
  };

  std::vector<std::unique_ptr<HTTPServer>> m_servers;
  std::vector<handler_reg> m_handlers;
  static http_server* g_server;
  static std::mutex g_mutex;
};

class http_request {
 public:
  connection& conn() { return m_connection; }
  const std::string& uri() const { return m_uri; }
  const std::string& method() const { return m_connection.m_req.getMethod(); }
  std::string peer() { return m_connection.m_req.clientAddress().host().toString(); }

  std::string uri_full() const { return m_connection.m_req.getURI(); }

  std::string get_header(const std::string& name) const;
  std::string get_header(const std::string& name, const std::string& default_value) const;
  void for_headers(const std::function<void(const std::string&, const std::string&)>& write) const;

  class variable_does_not_exist : public io_exception {
   public:
    variable_does_not_exist(const std::string& msg) : io_exception(msg) {}
  };

  std::string get_variable(const std::string& name);
  std::string get_variable(const std::string& name, const std::string& def);
  const std::map<std::string, std::string>& get_query_variables() { return m_query_variables; }
  std::string get_query() const;

  void send_continue();
  void send_status(int code, const std::string& message);
  // void send_status_line(const std::string& status_line);
  void send_header(const std::string& header, const std::string& value);
  void finish_headers();
  void finish_body(const std::string& body);
  const boost::smatch& uri_match() { return m_uri_match; }
  const boost::smatch& method_match() { return m_method_match; }

  bool is_ssl() const {
    // TODO(nils): make this more accurate.
    return false;
  }

  // If this request has credentials, return true and set the scheme, user & pw
  bool get_auth(std::string& scheme, std::string& user, std::string& password);

 private:
  friend struct http_server::handler_reg;

  http_request(HTTPServerRequest& req, HTTPServerResponse& resp);
  friend class RequestHandler;

  std::string find_variable(const std::string& name,
                            const std::function<std::string()>& no_query_string,
                            const std::function<std::string()>& var_not_found);

  connection m_connection;
  std::string m_uri;
  std::string m_query_str;
  boost::smatch m_uri_match;
  boost::smatch m_method_match;
  query_variables m_query_variables;
};
