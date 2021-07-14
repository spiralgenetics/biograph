#include "modules/pipeline/client_handler.h"
#include "modules/io/log.h" 
#include "modules/mapred/path.h" 
#include "modules/io/config.h" 

client_handler::client_handler(http_request& request)
	: rest_handler(request)
{}

void client_handler::get()
{
	std::string client_type{get_match_result(1)};
	SPLOG("client_handler::get> %s", client_type.c_str());

	path client_path(printstring("%s/bin/%s/spiral", CONF_CS(install_root), client_type.c_str()));
	if (!client_path.exists()) {
		m_request.send_status(404, "Client not found");
		write_output("Client could not be found.", "text/plain");
		return;
	}

	set_content_length(client_path.size());
	auto& output = set_output("application/binary", "spiral");
	io_copy(*client_path.read(), output);
}
