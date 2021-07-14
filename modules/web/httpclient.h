#pragma once

#include "modules/web/jsontypes.h"
#include "modules/io/io.h"

#include <string>
#include <map>
#include <list>
#include <memory>

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

typedef std::map<std::string, std::string> headers_type;

std::string make_client_url(
	const char* bind_list_var,
	const char* host_var,
	const char* port_var,
	const char* path = ""
);

class http_client
{
public:
	http_client(const std::string& base)
		: m_base(base)
	{}
	void set_cookie(const std::string& name, const std::string& value);
	void set_request_header(const std::string& key, const std::string& value);

	int do_request(const std::string& method, const std::string& url, const std::string& payload, std::string& result, bool spiral_request);

	int do_get(const std::string& url, std::string& result);
	int do_put(const std::string& url, const std::string& payload, std::string& result);
	int do_post(const std::string& url, const std::string& payload, std::string& result);
	int do_delete(const std::string& url);

	headers_type m_request_headers;
	headers_type m_response_headers;
	std::vector<HTTPCookie> m_cookies;

	size_t get_response_status_code() const { return m_last_status; }
	const std::string& get_response_status_message() const { return m_last_reason; }

private:
	void initialize_request_headers(); // set up essential request headers
	std::string m_base;
	size_t m_last_status = 0;
	std::string m_last_reason;
};

class json_client : public http_client
{
public:
	json_client(const std::string& base);
};

