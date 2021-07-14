#pragma once

#include "modules/web/restful.h" 
#include "modules/web/httpclient.h"

class forwarding_handler : public rest_handler
{
public:
	forwarding_handler(const std::string& endpoint, http_request& req);
	
	void get() override;
	void post() override;
	void put() override;
	void del() override;

private:
	void process_response(const std::string& response);
	bool rewrite_location_header(const std::string& location, std::string& new_location);

private:
	http_client m_http;
};

template <const std::string (*endpoint)()>
class forward_to : public forwarding_handler
{
public:
	forward_to(http_request& req) :
		forwarding_handler(endpoint(), req)
	{}
};
