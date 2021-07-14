#include "modules/web/forward_request_handler.h"
#include "modules/io/log.h"
#include <boost/regex.hpp>

forwarding_handler::forwarding_handler(const std::string& endpoint, http_request& req)
	: rest_handler(req)
	, m_http(endpoint)
{
	req.for_headers([&] (const std::string& name, const std::string& value) {
		m_http.set_request_header(name, value);
		// SPLOG_P(LOG_DEBUG, "forwarding_handler> %s: %s", name.c_str(), value.c_str());
	});
}

void forwarding_handler::get()
{
	// SPLOG_P(LOG_DEBUG, "forwarding_handler::get> %s", m_request.uri_full().c_str());
	std::string result;
	(void) m_http.do_get(m_request.uri_full(), result);
	process_response(result);
}

void forwarding_handler::post()
{
	// SPLOG_P(LOG_DEBUG, "forwarding_handler::post> %s", m_request.uri_full().c_str());
	std::string entity = read_entity(m_request, 64*1024*1024);
	// SPLOG_P(LOG_DEBUG, "POST entity:\n%s", entity.c_str());
	std::string result;
	(void) m_http.do_post(m_request.uri_full(), entity, result);
	process_response(result);
}

void forwarding_handler::put()
{
	// SPLOG_P(LOG_DEBUG, "forwarding_handler::put> %s", m_request.uri_full().c_str());
	std::string entity = read_entity(m_request, 64*1024*1024);
	// SPLOG_P(LOG_DEBUG, "PUT entity:\n%s", entity.c_str());
	std::string result;
	(void) m_http.do_put(m_request.uri_full(), entity, result);
	process_response(result);
}

void forwarding_handler::del()
{
	// SPLOG_P(LOG_DEBUG, "forwarding_handler::del> %s", m_request.uri_full().c_str());
	(void) m_http.do_delete(m_request.uri_full());
	process_response("");
}

void forwarding_handler::process_response(const std::string& result)
{
	// SPLOG_P(LOG_DEBUG, "forwarding_handler::process_response> entry");
	m_request.send_status(m_http.get_response_status_code(),
                          m_http.get_response_status_message());

	headers_type headers = m_http.m_response_headers;
	// SPLOG_P(LOG_DEBUG, "forwarding_handler::process_response> got %s headers", std::to_string(headers.size()).c_str());

	for(headers_type::iterator header = headers.begin(); header != headers.end(); header++) {

		// SPLOG_P(LOG_DEBUG, "forwarding_handler::process_response> %s: %s", header->first.c_str(), header->second.c_str());
		if (header->first == "Location") {
			std::string new_location;
			// SPLOG_P(LOG_DEBUG, "forwarding_handler::process_response> Response Location header: %s", header->second.c_str());
			if (rewrite_location_header(header->second, new_location)) {
				m_request.send_header(header->first, new_location);
			}
			else {
				m_request.send_header(header->first, header->second);
			}
		}
		else {
			m_request.send_header(header->first, header->second);
		}
	}
	if (headers.find("Content-Length") == headers.end()) {
		m_request.send_header("Content-Length", printstring("%lu", result.size()));
	}

	// cookies
	if(m_http.m_cookies.size() > 0) {
		std::vector<HTTPCookie>::iterator it = m_http.m_cookies.begin();
		for(; it != m_http.m_cookies.end(); ++it)
		{
			m_request.send_header("Set-Cookie", (*it).toString());
		}
	}

	m_request.finish_headers();

	// SPLOG_P(LOG_DEBUG, "forwarding_handler::process_response> Response body:\n%s", result.c_str());
	m_request.finish_body(result);
}

bool forwarding_handler::rewrite_location_header(
	const std::string& location,
	std::string& new_location)
{
	static boost::regex re_location("http://(?:[^/]+)(.*)");
	boost::smatch match;
	if (regex_match(location, match, re_location)) {
		new_location = printstring("%s", match[1].str().c_str());
		return true;
	}
	//SPLOG("Location did not match regex for Nginx");
	return false;
}
