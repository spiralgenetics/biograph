#include <assert.h>
#include <string.h>
#include <sys/prctl.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "Poco/Net/AcceptCertificateHandler.h"
#include "Poco/Net/Context.h"
#include "Poco/Net/HTTPBasicCredentials.h"
#include "Poco/Net/HTTPRequestHandler.h"
#include "Poco/Net/HTTPRequestHandlerFactory.h"
#include "Poco/Net/HTTPResponse.h"
#include "Poco/Net/HTTPServer.h"
#include "Poco/Net/HTTPServerRequest.h"
#include "Poco/Net/HTTPServerResponse.h"
#include "Poco/Net/NetSSL.h"
#include "Poco/Net/SSLManager.h"
#include "Poco/Net/SecureServerSocket.h"
#include "Poco/Net/ServerSocket.h"
#include "Poco/URI.h"
#include "Poco/Util/ServerApplication.h"

#include "base/base.h"
#include "modules/io/config.h"
#include "modules/io/hexdump.h"
#include "modules/io/log.h"
#include "modules/io/utils.h"
#include "modules/web/httpserver.h"

using namespace Poco::Net;
using namespace Poco::Util;
using namespace Poco;

std::mutex http_server::g_mutex;
http_server* http_server::g_server = nullptr;

class RequestHandler : public HTTPRequestHandler {
 public:
  RequestHandler(http_server* srv) : m_srv(srv) {}

  void handleRequest(HTTPServerRequest& request, HTTPServerResponse& response) override {
    http_request req(request, response);
    m_srv->handle(req);
  }

 private:
  http_server* m_srv;
};

class RequestHandlerFactory : public HTTPRequestHandlerFactory {
 public:
  RequestHandlerFactory(http_server* srv) : m_srv(srv) {}
  HTTPRequestHandler* createRequestHandler(const HTTPServerRequest& request) override {
    return new RequestHandler(m_srv);
  }

 private:
  http_server* m_srv;
};

int connection::base_read(char* buf, size_t len) {
  std::istream& istr = m_req.stream();
  auto ret = 0;
  // Read until len or EOF
  for (size_t i = 0; i < len; i++) {
    auto c = istr.get();
    if (c == EOF) {
      break;
    }
    buf[i] = c;
    ret++;
  }
  if (CONF(log_http_traffic)) {
    SPLOG_P(LOG_DEBUG, "connection::base_read> %d bytes", ret);
    if (ret > 0) {
      std::string chunk(buf, std::min(ret, 100));
      auto dump = hexdump(chunk);
      SPLOG_P(LOG_DEBUG, "\n%s", dump.c_str());
    }
  }
  return ret;
}

void connection::finish_headers() {
  CHECK(!m_send);
  m_send = &m_resp.send();
}

int connection::base_write(const char* buf, int len) {
  CHECK(m_send);

  m_send->write(buf, len);
  if (CONF(log_http_traffic)) {
    SPLOG_P(LOG_DEBUG, "connection::base_write> %d bytes", len);
    if (len > 0) {
      std::string chunk(buf, std::min(len, 100));
      auto dump = hexdump(chunk);
      SPLOG_P(LOG_DEBUG, "\n%s", dump.c_str());
    }
  }
  return len;
}

std::string http_request::get_header(const std::string& name) const {
  if (m_connection.m_req.has(name)) {
    return m_connection.m_req.get(name);
  }
  throw io_exception("No such header: " + name);
}

std::string http_request::get_header(const std::string& name, const std::string& def) const {
  return (m_connection.m_req.get(name, def));
}

void http_request::for_headers(
    const std::function<void(const std::string&, const std::string&)>& write) const {
  NameValueCollection& headers = m_connection.m_req;
  for (const auto& hdr : headers) {
    write(hdr.first, hdr.second);
  }
}

std::string http_request::find_variable(const std::string& name,
                                        const std::function<std::string()>& no_query_string,
                                        const std::function<std::string()>& var_not_found) {
  if (m_query_variables.empty()) {
    return no_query_string();
  }
  auto it = m_query_variables.find(name);
  if (it != m_query_variables.cend()) {
    return (*it).second;
  }
  return var_not_found();
}

std::string http_request::get_variable(const std::string& name) {
  return find_variable(
      name, []() -> std::string { throw variable_does_not_exist("No query string"); },
      [&name]() -> std::string {
        throw variable_does_not_exist(printstring("Failed to get variable %s", name.c_str()));
      });
}

std::string http_request::get_variable(const std::string& name, const std::string& def) {
  auto ret_def = [&def] { return def; };
  return find_variable(name, ret_def, ret_def);
}

std::string http_request::get_query() const { return m_query_str; }

void http_request::send_continue() { m_connection.m_resp.sendContinue(); }

void http_request::send_status(int code, const std::string& message) {
  m_connection.m_resp.setStatusAndReason(HTTPResponse::HTTPStatus(code), message);
}

void http_request::send_header(const std::string& header, const std::string& value) {
  // LOG("Response header: %s: %s", header.c_str(), value.c_str());
  m_connection.m_resp.set(header, value);
}

void http_request::finish_headers() { m_connection.finish_headers(); }

void http_request::finish_body(const std::string& body) {
  m_connection.write(body.c_str(), body.size());
}

http_request::http_request(HTTPServerRequest& req, HTTPServerResponse& resp)
    : m_connection(req, resp) {
  Poco::URI uri(req.getURI());
  m_uri = uri.getPath();
  m_query_str = uri.getQuery();
  for (const auto& param : uri.getQueryParameters()) {
    m_query_variables[param.first] = param.second;
  }

  // Strip out any double slashes in the URI.
  // TODO(nils): Figure out why we're generating double slashes in the first
  // place and remove this kludgery.
  for (;;) {
    auto pos = m_uri.find("//");
    if (pos == std::string::npos) {
      break;
    }
    m_uri = m_uri.substr(0, pos) + m_uri.substr(pos + 1);
  }
}

bool http_request::get_auth(std::string& scheme, std::string& user, std::string& password) {
  bool has_auth = m_connection.m_req.hasCredentials();
  std::string encoded;

  if (has_auth) {
    m_connection.m_req.getCredentials(scheme, encoded);

    HTTPBasicCredentials creds(m_connection.m_req);
    user = creds.getUsername();
    password = creds.getPassword();
  }

  return has_auth;
}

namespace {

static Context::Ptr get_ssl_context() {
  static bool ssl_initted = false;
  static Context::Ptr g_ssl_context;

  if (ssl_initted) {
    return g_ssl_context;
  }

  initializeSSL();

  SharedPtr<InvalidCertificateHandler> ptrCert;
  Context::Ptr ptrContext;

  g_ssl_context =
      new Context(Context::SERVER_USE, CONF_S(pem_file), CONF_S(ssl_certificates_chain),
                  "" /* no CA location specified; use default instead */, Context::VERIFY_NONE, 9,
                  true /* load default CAs */, "ALL:!ADH:!LOW:!EXP:!MD5:@STRENGTH");
  ptrCert = new AcceptCertificateHandler(false);
  SSLManager::instance().initializeServer(0, ptrCert, g_ssl_context);

  ssl_initted = true;
  return g_ssl_context;
}

}  // namespace

http_server& http_server::get() {
  std::lock_guard<std::mutex> l(g_mutex);
  if (g_server == nullptr) {
    g_server = new http_server();
  }
  return *g_server;
}

void http_server::register_handler(handler& handler, const std::string& uri_regex,
                                   const std::string& method_regex) {
  std::lock_guard<std::mutex> l(g_mutex);
  m_handlers.emplace_back(handler, uri_regex, method_regex);
}

void http_server::start(const bind_list_t& bl, const std::string& pem_path,
                        const std::string& ssl_certificates_chain_path) {
  try {
    std::lock_guard<std::mutex> l(g_mutex);
    CHECK(!bl.empty());
    for (size_t i = 0; i < bl.size(); i++) {
      auto port = bl[i].port;
      std::string ip = bl[i].ip;
      if (ip.empty()) {
        ip = "127.0.0.1";
      }
      ServerSocket sock;
      if (bl[i].ssl) {
        sock = SecureServerSocket(SocketAddress(ip, port), 64, get_ssl_context());
      } else {
        sock = ServerSocket(SocketAddress(ip, port));
      }
      auto new_http = new HTTPServer(new RequestHandlerFactory(this), sock, new HTTPServerParams);
      new_http->start();
      m_servers.emplace_back(new_http);
      SPLOG("Server listening on %s port %d%s", ip.c_str(), port, bl[i].ssl ? " (ssl)" : "");
    }
  } catch (const Poco::Exception& poco_except) {
    SPLOG("Got error starting server: %s", poco_except.message().c_str());
    throw std::runtime_error("Could not start server: " + poco_except.message());
  }
}

void http_server::stop() {
  std::lock_guard<std::mutex> l(g_mutex);
  SPLOG("Server shutting down");
  for (const auto& srv : m_servers) {
    srv->stop();
  }
  m_servers.clear();
}

void http_server::handle(http_request& request) {
  if (CONF(log_http_requests)) {
    SPLOG("%s -> %s %s", request.peer().c_str(), request.method().c_str(), request.uri().c_str());
  }
  for (auto& handler : m_handlers) {
    if (handler.match(request)) {
      request.conn().mark_valid();
      handler.m_handler.handle(request);
      return;
    }
  }
  SPLOG("HTTP error, code 404: Not Found -> %s", request.peer().c_str());
  request.send_status(404, "Not Found");
  request.send_header("Content-Type", "text/plain");
  request.finish_headers();
  request.conn().print(
      "Please use biograph to access this service. https://www.spiralgenetics.com/");
}

http_server::handler_reg::handler_reg(handler& handler, const std::string& uri_regex,
                                      const std::string& method_regex)
    : m_handler(handler) {
  try {
    uri = boost::regex(uri_regex, boost::regex_constants::extended);
    method = boost::regex(method_regex, boost::regex_constants::extended);
  } catch (const boost::regex_error& err) {
    throw io_exception("Invalid regex");
  }
}

bool http_server::handler_reg::match(http_request& request) {
  return boost::regex_match(request.uri(), request.m_uri_match, uri) &&
         boost::regex_match(request.m_connection.m_req.getMethod(), request.m_method_match, method);
}

struct error_response_json {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(success);
    FIELD(error);
    FIELD(status);
  }

  bool success = false;
  std::string error;
  unsigned status = 500;
};

void error_response(http_request& request, unsigned code, const std::string& msg) {
  SPLOG("HTTP error, code %d: %s -> %s", code, msg.c_str(), request.peer().c_str());
  request.send_status(code, msg);
  request.send_header("Content-Type", "application/json");
  error_response_json json;
  json.error = msg;
  json.status = code;
  auto body = json_serialize(json);
  request.send_header("Content-Length", printstring("%ld", body.size()));
  request.finish_headers();
  request.conn().write(body.data(), body.size());
}
