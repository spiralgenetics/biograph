#include "modules/web/httpsclient.h"
#include "modules/io/config.h"
#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/io/utils.h"

#include <boost/regex.hpp>
#include <deque>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>

#include "Poco/Exception.h"
#include "Poco/Net/AcceptCertificateHandler.h"
#include "Poco/Net/ConsoleCertificateHandler.h"
#include "Poco/Net/HTTPBasicCredentials.h"
#include "Poco/Net/HTTPCookie.h"
#include "Poco/Net/HTTPCredentials.h"
#include "Poco/Net/HTTPRequest.h"
#include "Poco/Net/HTTPResponse.h"
#include "Poco/Net/HTTPSClientSession.h"
#include "Poco/Net/NameValueCollection.h"
#include "Poco/Net/SSLManager.h"
#include "Poco/SharedPtr.h"
#include "Poco/StreamCopier.h"
#include "Poco/URI.h"

using Poco::StreamCopier;
using Poco::URI;
using Poco::Net::HTTPCookie;
using Poco::Net::HTTPMessage;
using Poco::Net::HTTPRequest;
using Poco::Net::HTTPResponse;
using Poco::Net::HTTPSClientSession;
using Poco::Net::NameValueCollection;

using namespace Poco::Net;
using namespace Poco;

class SSLInitializer {
 public:
  SSLInitializer() { Poco::Net::initializeSSL(); }

  ~SSLInitializer() { Poco::Net::uninitializeSSL(); }
};

class SpiralCertificateHandler : public InvalidCertificateHandler {
 public:
  SpiralCertificateHandler(bool handleErrorsOnServerSide)
      : InvalidCertificateHandler(handleErrorsOnServerSide){};

  void onInvalidCertificate(const void* pSender, VerificationErrorArgs& errorCert) override {
    SPLOG("SpiralCertificateHandler: cert failed verification");
    errorCert.setIgnoreError(false);
  }
};

int https_client::do_request(const std::string& method, const std::string& url,
                             const std::string& payload, std::string& result) {
  if (CONF(log_http_traffic)) {
    SPLOG_P(LOG_DEBUG, "https_client::do_request> %s %s", method.c_str(), (m_base + url).c_str());
    SPLOG_P(LOG_DEBUG, "https_client::do_request> %s", payload.c_str());
  }
  m_last_status = 520;  // Unknown Error

  // Disable keepalive
  m_request_headers["Connection"] = "close";

  try {
    URI uri(m_base + url);

    SSLInitializer sslInitializer;

    SharedPtr<InvalidCertificateHandler> ptrCert;
    Context::Ptr ptrContext;
    std::string cipher_list = "ALL:!ADH:!LOW:!EXP:!MD5:@STRENGTH";

    if (m_skip_validation) {
      // https://pocoproject.org/docs/Poco.Net.Context.html
      // Context(usage, privateKeyFile, certificateFile, caLocation, VerificationMode,
      // verificationDepth, loadDefaultCAs, cipherList)
      ptrContext =
          new Context(Context::CLIENT_USE, "", "", "", Context::VERIFY_NONE, 9, false, cipher_list);
      ptrCert = new AcceptCertificateHandler(false);
    } else {
      // Only load the default CAs if we weren't given one in m_ca.
      ptrContext = new Context(Context::CLIENT_USE, "", "", "", Context::VERIFY_STRICT, 9,
                               m_ca.empty(), cipher_list);
      if (not m_ca.empty()) {
        std::stringstream ca_ss;
        ca_ss << m_ca;
        ptrContext->addCertificateAuthority(X509Certificate(ca_ss));
      }
      // Disable hostname check, https://pocoproject.org/docs/Poco.Net.Context.html#26288
      ptrContext->enableExtendedCertificateVerification(false);
      ptrCert = new SpiralCertificateHandler(false);
    }
    SSLManager::instance().initializeClient(0, ptrCert, ptrContext);

    HTTPSClientSession session(uri.getHost(), uri.getPort());
    HTTPRequest request(method, uri.getPathAndQuery(), HTTPMessage::HTTP_1_1);
    HTTPResponse response;

    // Auth
    HTTPBasicCredentials creds(m_user, m_password);
    creds.authenticate(request);

    // Set headers
    for (headers_type::iterator it = m_request_headers.begin(); it != m_request_headers.end();
         it++) {
      request.set(it->first, it->second);
    }

    // Set cookies (if any)
    if (m_cookies.size() > 0) {
      NameValueCollection nvc;
      std::vector<HTTPCookie>::iterator it = m_cookies.begin();
      for (; it != m_cookies.end(); ++it) {
        nvc.add((*it).getName(), (*it).getValue());
      }
      request.setCookies(nvc);
    }

    // Poco requires chunked encoding or manually setting the Content-Length header.
    request.setContentLength(payload.length());

    // One minute timeout.
    session.setTimeout(Poco::Timespan(60, 0));

    // Force keepalive to false
    session.setKeepAlive(false);

    session.sendRequest(request) << payload;
    std::istream& rs = session.receiveResponse(response);
    StreamCopier::copyToString(rs, result);

    m_last_status = response.getStatus();
    m_last_reason = response.getReason();

    // if (response.getStatus() != HTTPResponse::HTTP_UNAUTHORIZED) ...

    // Save response headers
    m_response_headers.clear();
    NameValueCollection::ConstIterator it = response.begin();
    for (; it != response.end(); ++it) {
      m_response_headers[it->first] = it->second;
    }

    // Keep cookies (if any)
    response.getCookies(m_cookies);
  } catch (std::exception& e) {
    SPLOG_P(LOG_DEBUG, "https_client::do_request> exception: %s", e.what());
    throw(std::runtime_error(e.what()));
  }

  if (CONF(log_http_traffic)) {
    SPLOG_P(LOG_DEBUG, "https_client::do_request> status: %lu", m_last_status);
    SPLOG_P(LOG_DEBUG, "https_client::do_request> reason: '%s'", m_last_reason.c_str());
    SPLOG_P(LOG_DEBUG, "https_client::do_request> result: '%s'", result.c_str());
  }

  return m_last_status;
}

int https_client::do_get(const std::string& url, std::string& result) {
  std::string payload;
  return do_request(HTTPRequest::HTTP_GET, url, payload, result);
}

int https_client::do_put(const std::string& url, const std::string& payload, std::string& result) {
  return do_request(HTTPRequest::HTTP_PUT, url, payload, result);
}

int https_client::do_post(const std::string& url, const std::string& payload, std::string& result) {
  return do_request(HTTPRequest::HTTP_POST, url, payload, result);
}

int https_client::do_delete(const std::string& url) {
  std::string result;
  std::string payload;
  return do_request(HTTPRequest::HTTP_DELETE, url, payload, result);
}

void https_client::set_credentials(const std::string& user, const std::string& password) {
  m_user = user;
  m_password = password;
}
