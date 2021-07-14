#include "modules/mapred/resource_manager.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/name_generator.h"
#include "modules/io/config.h"
#include "modules/io/uuid.h"

#include <boost/filesystem.hpp>

#include <sys/file.h>
#include <sys/fcntl.h>
#include <sys/stat.h>

class mmap_readable : public reset_readable
{
public:
	mmap_readable(const char* buf, size_t size)
		: m_buf(buf)
		, m_size(size)
	{}

	void reset() override
	{
		m_offset = 0;
	}

	size_t read(char* buf, size_t len) override
	{
		len = std::min(len, m_size - m_offset);
		std::copy_n(m_buf + m_offset, len, buf);
		m_offset += len;
		return len;
	}

private:
	const char* m_buf;
	size_t m_size;
	size_t m_offset = 0;
};

resource_manager::resource_manager()
{
	path storage_root(CONF_S(storage_root));
	if (storage_root.type() == path::FILE) {
		m_direct = true;
	}
}

resource_manager::resource_manager(bool direct)
	: m_direct(direct)
{}

void resource_manager::create_resource(mmap_buffer& out, size_t size)
{
	path root(CONF_S(resources_root));
	root.mkdir();
	std::string dir = root.bare_path();
	auto id = make_uuid();
	std::string path = printstring("%s/%s", dir.c_str(), id.c_str());

	SPLOG_P(LOG_DEBUG, "resource_manager::create_resource> size = %lu, %s", size, path.c_str());

	if (!m_direct) {
		make_space(size);
	}

	out.open(path, size);
	out.set_uuid(id);
}

void resource_manager::write_resource(
	manifest& out,
	mmap_buffer& in,
	const path& root,
	const std::string& prefix,
	const progress_handler_t& progress)
{
	if (m_direct) {
		path root(CONF_S(path_bulkdata));
		root.mkdir();

		path src(in.path());
		path dest(root.append(in.get_uuid()));

		SPLOG_P(LOG_DEBUG, "resource_manager::write_resource> %s -> %s",
			src.bare_path().c_str(),
			dest.bare_path().c_str()
		);

		path::move(src, dest);

		file_info fi;
		fi.file = dest;
		fi.size = dest.size();
		fi.num_records = 0;
		out.add(fi, 0);
		progress(1.0);

		return;
	}

	SPLOG_P(LOG_DEBUG, "resource_manager::write_resource> %s", in.path().c_str());
	std::string pname = in.path();

	const size_t chunk_size = 64 * 1024 * 1024;
	int chunk_num = 0;
	for (size_t i = 0; i < in.size(); i += chunk_size) {
		progress(double(i) / double(in.size()));
		size_t cur_size = std::min(chunk_size, in.size() - i);
		mmap_readable mr(in.buffer() + i, cur_size);
		auto name = printstring("%s_%d", prefix.c_str(), chunk_num++);
		auto tmppath = root.append(printstring("%03ld_%s_%s",
			random() % 1000,
			make_uuid().c_str(),
			name.c_str())
		);
		tmppath.write_inverted(mr, cur_size)->wait();
		file_info fi;
		fi.file = tmppath;
		fi.size = cur_size;
		fi.num_records = 0;
		out.add(fi, 0);
	}
	out.metadata().set(meta::ns::internal, "resource_uuid", in.get_uuid());
	in.close();
	chmod(pname.c_str(), 0444);
	in.open(pname, mmap_buffer::mode::copy_on_write);
	progress(1.0);
}

static
void write_complete(int fd, const char* buf, size_t len)
{
	size_t offset = 0;
	while (offset < len) {
		ssize_t r = write(fd, buf + offset, len - offset);
		if (r < 0 && errno == EINTR) {
			continue;
		}
		if (r < 0) {
			throw io_exception("Error in write_complete");
		}
		offset += r;
	}
}

void resource_manager::read_resource(
	mmap_buffer& out,
	const manifest& in,
	const progress_handler_t& progress)
{
	if (m_direct) {
		auto first = in.begin();
		SPLOG_P(LOG_DEBUG, "resource_manager::read_resource> %s", first->file.bare_path().c_str());
		out.open(first->file.bare_path(), mmap_buffer::mode::copy_on_write);
		progress(1.0);
		return;
	}

	make_space(in.get_size());
	std::string uuid = in.metadata().get<std::string>(meta::ns::internal, "resource_uuid");
	auto dir = CONF_S(resources_root);
	std::string path = printstring("%s/%s", dir.c_str(), uuid.c_str());
	SPLOG_P(LOG_DEBUG, "resource_manager::read_resource> %s", path.c_str());
	// Make make/download resource
	while (true) {
		//SPLOG("Trying open");
		// TODO: Use a better file locking mechanism here. This scheme fails if the process
		// runs as root, since it can always open the file for writing.
		int fd = open(path.c_str(), O_RDWR | O_CREAT, 0664);
		if (fd < 0 && errno == EACCES) {
			SPLOG("resource_manager::read_resource> Open in write mode failed %s", path.c_str());
			break; // Probably alrady RO, give the mapping a try
		}
		if (fd < 0) {
			throw io_exception("resource_manager: Unable to create file"); // Otherwise, bad news
		}
		int r = flock(fd, LOCK_EX | LOCK_NB);
		if (r < 0 && errno == EWOULDBLOCK) {
			// Some else is downloading, try again
			SPLOG("resource_manager::read_resource> flock() failed %s", path.c_str());
			close(fd);
			progress(0.0);
			sleep(3);
			continue;
		}
		if (r < 0) {
			throw io_exception("resource_manager: Unable to lock file");
		}
		// I've got the lock, do the download
		SPLOG("resource_manager::read_resource> Loading data from manifest");
		manifest_reader mr(in);
		const size_t k_bufsize = 16*1024;
		char buf[k_bufsize];
		size_t tot_size = in.get_size();
		size_t cur_size = 0;
		while (true) {
			size_t r = mr.read(buf, k_bufsize);
			if (r == 0) {
				break;
			}
			write_complete(fd, buf, r);
			cur_size += r;
			progress(double(cur_size) / double(tot_size));
		}
		//SPLOG("Loaded the data, truncating\n");
		r = ftruncate(fd, tot_size);
		if (r < 0) {
			throw io_exception("resource_manager: Unable to truncate file");
		}
		// Mark read only (download done), and close
		//SPLOG("Changing mode");
		fchmod(fd, 0444);
		close(fd);
	}
	SPLOG_P(LOG_DEBUG, "resource_manager::read_resource> Doing RO open");
	out.open(path, mmap_buffer::mode::copy_on_write);
	progress(1.0);
}

size_t resource_manager::free_space()
{
	namespace fs = boost::filesystem;
	fs::path dir = CONF_S(resources_root);
	auto si = fs::space(dir);
	size_t slop = CONF_T(size_t, resource_quota_slop);
	if (si.available < slop) {
		return 0;
	}
	return si.available - slop;
}

void resource_manager::make_space(size_t space)
{
	namespace fs = boost::filesystem;
	try {
		fs::path dir = CONF_S(resources_root);
		fs::create_directory(dir);
		while (free_space() < space) {
			fs::path which;
			std::time_t oldest = std::numeric_limits<std::time_t>::max();
			for (auto it = fs::directory_iterator(dir); it != fs::directory_iterator(); ++it) {
				std::time_t last_write = fs::last_write_time(*it);
				if (last_write < oldest) {
					oldest = last_write;
					which = *it;
				}
			}
			if (which.empty()) {
				throw io_exception("No room for resource!");
			}
			fs::remove(which);
		}
	} catch (const fs::filesystem_error& err) {
		throw io_exception(std::string("Error during make_space: ") + err.what());
	}
}
