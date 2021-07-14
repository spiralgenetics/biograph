#include "modules/pipeline/version.h"
#include "modules/io/json_transfer.h"
#include "modules/io/version.h"

version_handler::version_handler(http_request& req)
	: easy_rest_handler(req)
{
}

std::string version_handler::easy_get()
{
	return "unimplemented";
}
