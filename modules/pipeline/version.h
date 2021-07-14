#pragma once

#include "modules/web/restful.h"

class version_handler : public easy_rest_handler
{
public:
	version_handler(http_request& req);
	std::string easy_get() override;
};
