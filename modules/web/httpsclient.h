#pragma once

#include "modules/io/io.h"
#include "modules/web/httpclient.h"
#include "modules/web/jsontypes.h"

#include <list>
#include <map>
#include <memory>
#include <string>

#include "Poco/Net/HTTPCookie.h"
#include "Poco/Net/HTTPCredentials.h"
#include "Poco/Net/HTTPRequest.h"
#include "Poco/Net/HTTPResponse.h"
#include "Poco/Net/HTTPSClientSession.h"
#include "Poco/Net/NameValueCollection.h"
#include "Poco/StreamCopier.h"
#include "Poco/URI.h"

using Poco::StreamCopier;
using Poco::URI;
using Poco::Net::HTTPCookie;
using Poco::Net::HTTPCredentials;
using Poco::Net::HTTPMessage;
using Poco::Net::HTTPRequest;
using Poco::Net::HTTPResponse;
using Poco::Net::HTTPSClientSession;
using Poco::Net::NameValueCollection;

class https_client : public http_client {
 public:
  https_client(const std::string& base) : http_client(base), m_base(base) {}

  // Constructor for including a trusted certificate authority string
  https_client(const std::string& base, const std::string& ca) : http_client(base), m_base(base), m_ca(ca) {}

  int do_request(const std::string& method, const std::string& url, const std::string& payload,
                 std::string& result);

  int do_get(const std::string& url, std::string& result);
  int do_put(const std::string& url, const std::string& payload, std::string& result);
  int do_post(const std::string& url, const std::string& payload, std::string& result);
  int do_delete(const std::string& url);

  bool skip_validation(bool flag) { m_skip_validation = flag; return m_skip_validation; }
  bool skip_validation() { return m_skip_validation; }

  void set_credentials(const std::string& user, const std::string& password);

 private:
  std::string m_base;
  std::string m_ca;
  size_t m_last_status = 0;
  std::string m_last_reason;
  bool m_skip_validation = false;
  std::string m_user;
  std::string m_password;
};
