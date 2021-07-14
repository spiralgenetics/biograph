#include "modules/io/utils.h"
#include "modules/io/io.h"
#include "modules/io/config.h"
#include "modules/mapred/path.h"
#include "modules/mapred/temp_file.h"

#include <vector>
#include <cstring>
#include <unistd.h>

#include <boost/filesystem.hpp>

scoped_temp_file::scoped_temp_file()
	: scoped_temp_file(::path(CONF_S(temp_root) + "/spiral-XXXXXX").bare_path())
{
}

scoped_temp_file::scoped_temp_file(const std::string& tmpl)
{
	boost::filesystem::path temp_path(tmpl);
	if (! boost::filesystem::is_directory(temp_path.parent_path())) {
		boost::filesystem::create_directories(temp_path.parent_path());
	}

	std::vector<char> buf(tmpl.begin(), tmpl.end());
	buf.push_back('\0');
	m_fd = ::mkstemp(buf.data());
	if (m_fd == -1) {
		throw io_exception(printstring("scoped_temp_file> ::mkstemp(%s) failed: %s",
			buf.data(), strerror(errno)));
	}
	m_path = buf.data();
}

scoped_temp_file::~scoped_temp_file()
{
	::close(m_fd);
	if (!m_path.empty()) {
		::unlink(m_path.c_str());
	}
}
