#pragma once
#include "modules/web/restful.h"

class client_handler : public rest_handler
{
public:
	client_handler(http_request& req);
	void get() override;
};
