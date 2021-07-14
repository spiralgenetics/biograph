#include "modules/web/restful.h"
#include "modules/io/config.h"
#include "modules/io/log.h"
#include "modules/io/utils.h"
#include "modules/web/chunked_encoding.h"
#include "modules/web/jsontypes.h"

#include <thread>

class master_rest_handler : public handler {
 public:
  master_rest_handler(create_handler_t create_handler) : m_create_handler(create_handler) {}

  void handle(http_request& request) override;

 private:
  create_handler_t m_create_handler;
};

void master_rest_handler::handle(http_request& request) {
  try {
    std::string method = request.method();
    std::unique_ptr<rest_handler> handler(m_create_handler(request));

    if (!handler->auth()) {
      throw forbidden("Not authorized");
    }

    if (method == "GET") {
      handler->get();
    } else if (method == "PUT") {
      handler->put();
    } else if (method == "POST") {
      handler->post();
    } else if (method == "DELETE") {
      handler->del();
    } else if (method == "PATCH") {
      handler->patch();
    } else if (method == "OPTIONS") {
      handler->options();
    } else {
      throw method_not_allowed(method);
    }
  } catch (const rest_exception& re) {
    error_response(request, re.get_errcode(), re.message());
  } catch (const io_exception& io) {
    SPLOG("master_rest_handler::handle> io_exception: %s", io.message().c_str());
    error_response(request, 500, "Internal server error: " + io.message());
  }
}

rest_handler::rest_handler(http_request& request)
    : m_request(request), m_chunker(m_request.conn()) {}

std::string rest_handler::get_cookie(const std::string& name) const {
  std::string all_cookies;
  try {
    all_cookies = " " + m_request.get_header("Cookie");
  } catch (const io_exception& io) {
    return "";
  }

  // TODO:  This does not properly deal with qouted cookies
  std::string key = " " + name + "=";
  std::string::size_type r = all_cookies.find(key);
  if (r == std::string::npos) {
    return "";
  }
  std::string::size_type value_start = r + key.size();
  std::string::size_type value_end = all_cookies.find(";", value_start);
  if (value_start != std::string::npos) {
    return all_cookies.substr(value_start, value_end - value_start);
  }
  return all_cookies.substr(value_start);
}

void rest_handler::set_status_code(int status, const std::string& message) {
  m_request.send_status(status, message);
}

void rest_handler::set_cookie(const std::string& cookie, const std::string& value, int max_age) {
  std::string cookiestr =
      printstring("%s=%s; path=/; max-age=%d", cookie.c_str(), value.c_str(), max_age);
  m_request.send_header("Set-Cookie", cookiestr);
}

void rest_handler::set_content_length(size_t length) {
  std::string lenstr = printstring("%lu", (unsigned long)length);
  m_request.send_header("Content-Length", lenstr);
}

writable& rest_handler::set_output(const std::string& content_type, const std::string& filename,
                                   bool use_chunked_transfer_encoding) {
  m_request.send_header("Content-Type", content_type);
  if (filename != "") {
    m_request.send_header("Content-Disposition", "attachment; filename=" + filename);
  }
  if (use_chunked_transfer_encoding) {
    m_request.send_header("Transfer-Encoding", "chunked");
  }
  if (content_type.find(JSONTYPE) != std::string::npos) {
    m_request.send_header("Cache-Control", "no-cache, no-store, must-revalidate");
    m_request.send_header("Pragma", "no-cache");
    m_request.send_header("Expires", "0");
  }
  m_request.finish_headers();
  if (use_chunked_transfer_encoding) {
    return m_chunker;
  }
  return m_request.conn();
}

void rest_handler::write_output(const std::string& content, const std::string& content_type) {
  set_content_length(content.size());
  writable& out = set_output(content_type);
  out.write(content.c_str(), content.size());
}

void easy_rest_handler::get() {
  std::string response = easy_get();
  if (response.empty()) {
    set_status_code(204);
  }
  write_output(response, content_type());
}

void easy_rest_handler::put() {
  std::string input = read_entity(m_request, max_size());
  if (easy_put(input)) {
    set_status_code(201);
  }
  write_output("", "text/plain");
}

void easy_rest_handler::post() {
  std::string input = read_entity(m_request, max_size());
  std::string output = easy_post(input);
  write_output(output, content_type());
}

void easy_rest_handler::del() {
  bool response = easy_del();
  write_output(response ? "true" : "false", JSONTYPE);
}

void easy_rest_handler::patch() {
  std::string input = read_entity(m_request, max_size());
  std::string output = easy_patch(input);
  write_output(output, content_type());
}

void easy_rest_handler::options() {
  std::string input = read_entity(m_request, max_size());
  std::string output = easy_options(input);
  write_output(output, content_type());
}

void register_handler(const std::string& path, create_handler_t create_handler) {
  master_rest_handler* h = new master_rest_handler(create_handler);
  http_server::get().register_handler(*h, path, "(.*)");
}

void cat_nap() { std::this_thread::sleep_for(std::chrono::milliseconds(500)); }

void run_restful_server(const bind_list_t& bind_list, const std::string& pem_path,
                        const std::string& ssl_certificates_chain_path,
                        const std::string& fork_mode, bool block) {
  // Set config to match supplied pem and cert
  Config::set("pem_file", pem_path);
  Config::set("ssl_certificates_chain", ssl_certificates_chain_path);

  if (fork_mode == "thread") {
    http_server::get().start(bind_list, pem_path, ssl_certificates_chain_path);
    while (block) {
      cat_nap();
    }
    return;
  } else if (fork_mode == "fork") {
    pid_t pid = fork();
    if (pid != 0) {
      while (block) {
        cat_nap();
      }
      return;
    }
    // double fork
    pid = fork();
    if (pid != 0) {
      // child exits
      exit(0);
    }
    // grandchild is adopted by init (to avoid zombies)
    http_server::get().start(bind_list, pem_path, ssl_certificates_chain_path);
    while (true) {
      cat_nap();
    }
  } else {
    throw std::runtime_error("Invalid fork_mode: " + fork_mode);
  }
}

std::string read_entity(http_request& request, size_t max_size) {
  std::string slen;
  try {
    slen = request.get_header("Content-Length");
  } catch (const io_exception& io) {
    throw rest_exception("Length required", 411);
  }
  size_t len = atol(slen.c_str());
  if (len > max_size) {
    throw rest_exception("Request too large", 413);
  }

  std::string result;
  result.resize(len);
  size_t rlen = request.conn().read(&result[0], len);
  if (rlen != len) {
    throw rest_exception("Failed to read content", 400);
  }

  return result;
}
