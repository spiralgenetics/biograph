#include "modules/web/httpclient.h"
#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/io/utils.h"
#include "modules/io/make_unique.h"
#include "modules/io/config.h"
#include "modules/web/httpserver.h"

#include <vector>
#include <fstream>
#include <deque>
#include <utility>
#include <boost/regex.hpp>

#include "Poco/Net/HTTPClientSession.h"
#include "Poco/Net/HTTPRequest.h"
#include "Poco/Net/HTTPResponse.h"
#include "Poco/Net/HTTPCookie.h"
#include "Poco/Net/NameValueCollection.h"
#include "Poco/URI.h"
#include "Poco/StreamCopier.h"

using Poco::Net::HTTPClientSession;
using Poco::Net::HTTPRequest;
using Poco::Net::HTTPResponse;
using Poco::Net::HTTPMessage;
using Poco::Net::HTTPCookie;
using Poco::Net::NameValueCollection;
using Poco::StreamCopier;
using Poco::URI;

void http_client::set_cookie(const std::string& name, const std::string& value)
{
	m_cookies.push_back(HTTPCookie(name, value));
}

void http_client::set_request_header(const std::string& key, const std::string& value)
{
	m_request_headers[key] = value;
	return;
}

int http_client::do_request(const std::string& method, const std::string& url, const std::string& payload, std::string& result, bool spiral_request = true) {
	m_last_status = 520; // Unknown Error

	if(spiral_request) {
		// Disable keepalive
		m_request_headers["Connection"] = "close";
	}

	// SPLOG_P(LOG_DEBUG, "http_client::do_request> %s %s", method.c_str(), (m_base + url).c_str());

	try {
		URI uri(m_base + url);

		HTTPClientSession session(uri.getHost(), uri.getPort());
		HTTPRequest request(method, uri.getPathAndQuery(), HTTPMessage::HTTP_1_1);
		HTTPResponse response;

		// Set headers
		for(headers_type::iterator it = m_request_headers.begin();
			it != m_request_headers.end(); it++) {
			request.set(it->first, it->second);
		}

		// Set cookies (if any)
		if(m_cookies.size() > 0) {
			NameValueCollection nvc;
			std::vector<HTTPCookie>::iterator it = m_cookies.begin();
			for(; it != m_cookies.end(); ++it)
			{
				nvc.add((*it).getName(), (*it).getValue());
			}
			request.setCookies(nvc);
		}

		// Poco requires chunked encoding or manually setting the Content-Length header.
		request.setContentLength(payload.length());

		// One minute timeout.
		session.setTimeout(Poco::Timespan(60,0));

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
		for(; it != response.end(); ++it) {
			m_response_headers[it->first] = it->second;
		}

		// Keep cookies (if any)
		response.getCookies(m_cookies);
	}
	catch (std::exception &e) {
		SPLOG_P(LOG_DEBUG, "http_client::do_request> exception: %s", e.what());
		throw(std::runtime_error(e.what()));
	}

	if(CONF(log_http_traffic)) {
		SPLOG_P(LOG_DEBUG, "http_client::do_request> status: %lu", m_last_status);
		SPLOG_P(LOG_DEBUG, "http_client::do_request> reason: '%s'", m_last_reason.c_str());
		SPLOG_P(LOG_DEBUG, "http_client::do_request> result: '%s'", result.c_str());
	}

	return m_last_status;
}

int http_client::do_get(const std::string& url, std::string& result)
{
	std::string payload;
	return do_request(HTTPRequest::HTTP_GET, url, payload, result);
}

int http_client::do_put(const std::string& url, const std::string& payload, std::string& result)
{
	return do_request(HTTPRequest::HTTP_PUT, url, payload, result);
}

int http_client::do_post(const std::string& url, const std::string& payload, std::string& result)
{
	return do_request(HTTPRequest::HTTP_POST, url, payload, result);
}

int http_client::do_delete(const std::string& url)
{
	std::string result;
	std::string payload;
	return do_request(HTTPRequest::HTTP_DELETE, url, payload, result);
}

std::string make_client_url(
	const char* bind_list_var,
	const char* host_var,
	const char* port_var,
	const char* path)
{
	bind_list_t bind_list{Config::instance().get<bind_list_t>(bind_list_var)};
	if (bind_list.empty()) {
			throw io_exception(printstring("Missing bind_list named: %s", bind_list_var));
	}

	if (bind_list[0].ip.empty()) {
			// default to 127.0.0.1 if the bind_list doesn't specify an IP
			bind_list[0].ip = "127.0.0.1";
	}

	auto host = getenv_str(host_var, bind_list[0].ip);
	auto port = getenv_int(port_var, bind_list[0].port);
	return printstring("http://%s:%d%s", host.c_str(), port, path);
}

json_client::json_client(const std::string& base)
	: http_client(base)
{
	m_request_headers["Content-Type"] = JSONTYPE;
	m_request_headers["Expect"] = "";
}


