#include <fcntl.h>

#include <boost/format.hpp>

#include "modules/io/io.h"
#include "modules/io/file_wrapper.h"


file_wrapper::file_wrapper(const std::string& path, int open_flags, mode_t mode)
{
	open(path, open_flags, mode);
}


void file_wrapper::open(const std::string& path, int open_flags, mode_t mode)
{
	if (is_open()) {
		throw io_exception(boost::format("Trying to open already-open file %1%") % path);
	}

	m_fd = ::open(path.c_str(), open_flags, mode);
	if (m_fd < 0) {
		if ((errno == 17) && (open_flags & O_EXCL)) {
			throw io_exception(boost::format(
				"file_wrapper::open> open() failed for path %1% with flags %2$#x, errno = %3% "
				"attempting to open a file exclusively (O_EXCL) when it already exists.")
				% path % open_flags % errno
			);
		} else {
			throw io_exception(boost::format(
				"file_wrapper::open> open() failed for path %1% with flags %2$#x, errno = %3%")
				% path % open_flags % errno
			);
		}
	}
	
	m_path = path;
}


void file_wrapper::close() noexcept
{
	if (m_fd == -1) {
		return;
	}

	::close(m_fd);
	m_fd = -1;
	m_path.clear();
}
