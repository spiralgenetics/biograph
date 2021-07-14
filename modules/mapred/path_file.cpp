#include "modules/mapred/path_impl.h"
#include "modules/io/file_io.h"
#include "modules/io/hash_io.h"
#include "modules/io/zip.h"
#include "modules/io/make_unique.h"
#include "modules/io/base64.h"

#define BOOST_FILESYSTEM_VERSION 3
#include <boost/version.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <sys/stat.h>
#include <future>

static const boost::regex file_regex("file://(/.*)");

static void make_subdir(const std::string& str)
{
	boost::filesystem::path p(str);
	try {
		boost::filesystem::create_directories(p.branch_path());
	}
	catch (const std::exception & ex) {
		// Ignore
	}
}

class file_write_waiter : public waiter
{
public:
	file_write_waiter(reset_readable& source, const std::string& path)
		: m_source(source)
		, m_path(path)
		, m_future(std::async(std::launch::async, std::bind(&file_write_waiter::run, this)))
	{}

	~file_write_waiter()
	{
		wait();
	}

	std::string wait() override
	{
		if (m_future.valid()) {
			m_md5_base64 = m_future.get();
		}
		return m_md5_base64;
	}

private:
	std::string run()
	{
		make_subdir(m_path);
		file_writer fw(m_path);
		md5_hash_writer hw;
		multi_writer writer(&fw, &hw);
		io_copy(m_source, writer);
		hw.finish();
		auto digest = hw.digest();
		return base64_encode(digest.data(), digest.size());
	}

private:
	reset_readable& m_source;
	std::string m_path;
	std::future<std::string> m_future;
	std::string m_md5_base64;
};

path_file_impl::path_file_impl(const std::string& url)
{
	boost::smatch parts;
	if (boost::regex_match(url, parts, file_regex)) {
		m_path = parts[1];		
	}
	else {
		m_path = url;
	}
}

path_impl* path_file_impl::clone() const
{
	return new path_file_impl(url());
}

std::string path_file_impl::url() const
{
	if (m_path[0] == '/') {
		return std::string("file://") + m_path;
	}
	return m_path;
}

std::unique_ptr<readable> path_file_impl::read() const
{
	return make_unique<file_reader>(m_path);
}

std::unique_ptr<writable> path_file_impl::write(const path_write_options& options) const
{
	make_subdir(m_path);
	return make_unique<file_writer>(m_path);
}

std::unique_ptr<waiter> path_file_impl::write_inverted(reset_readable& source, size_t size, const path_write_options& options) const
{
	return make_unique<file_write_waiter>(source, m_path);
}

void path_file_impl::move(const path& src, const path& dest) const
{
	std::rename(m_path.c_str(), dest.bare_path().c_str());	
}

void path_file_impl::copy(const path& src, const path& dest, const path_write_options& options) const
{
	auto reader = read();
	auto writer = dest.write();
	io_copy(*reader, *writer);
}

path::exist_enum path_file_impl::exists() const 
{
	if (boost::filesystem::is_regular_file(m_path)) {
		return path::e_file;
	}
	else if (boost::filesystem::is_directory(m_path)) {
		return path::e_directory;
	}
	return path::e_no_exist;
}

std::time_t path_file_impl::modify_time() const
{
	try {
		return boost::filesystem::last_write_time(m_path);
	}
	catch (...) {
		throw io_exception("Couldn't find time for " + m_path);
	}
}

size_t path_file_impl::size() const
{
	try {
		return boost::filesystem::file_size(m_path);
	}
	catch (...) {
		throw io_exception("Couldn't get size for " + m_path);
	}
}

std::vector<std::string> path_file_impl::list() const
{
	std::vector<std::string> out;
	if (!boost::filesystem::is_directory(m_path)) {
		throw io_exception(printstring("Trying to list a non-directory: %s", m_path.c_str()));
	}

	typedef boost::filesystem::directory_iterator iterator;
	for (iterator it(m_path); it != iterator(); ++it) {
#if BOOST_VERSION <= 104300
		out.push_back(it->path().filename());
#else
		out.push_back(it->path().filename().string());
#endif
	}
	return out;
}

void path_file_impl::mkdir() const
{
	boost::filesystem::create_directories(m_path);
}

bool path_file_impl::rm() const
{
	int rc = unlink(m_path.c_str());
	if (rc == 0) {
		return true;
	}
	if (rc == -1 && errno == ENOENT) {
		return false;
	}
	throw io_exception(printstring("Unable to delete file: %s\n", strerror(errno)));
}

bool path_file_impl::rmdir() const
{
	return boost::filesystem::remove(m_path);
}

static 
void walk_rec(const boost::filesystem::path& p, const path::walker_f& fn) 
{
	struct stat sb;
	if (stat(p.string().c_str(), &sb) < 0) {
		throw io_exception(printstring("Unable to stat file in path_file_impl::walk: %s", 
			p.string().c_str()
		));
	}
	path ps(p.string());
	if ((sb.st_mode & S_IFMT) == S_IFDIR) {
		fn(path::walk_params(path::w_dir_enter, ps, sb.st_mtime, sb.st_size));
		boost::filesystem::directory_iterator it(p), itEnd;
		for (; it != itEnd; ++it) {
			walk_rec(it->path(), fn);
		}
		fn(path::walk_params(path::w_dir_leave, ps, sb.st_mtime, sb.st_size));
	} else {
		fn(path::walk_params(path::w_file, ps, sb.st_mtime, sb.st_size));
	}
}

void path_file_impl::walk(path::walker_f fn) const
{
	if (exists() == path::e_no_exist) {
		return;
	}
	walk_rec(boost::filesystem::path(m_path), fn);
}
