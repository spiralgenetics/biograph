#include <boost/algorithm/string/predicate.hpp>
#include <stack>

#include "modules/io/config.h"
#include "modules/io/log.h"
#include "modules/mapred/path.h"
#include "modules/mapred/path_impl.h"
#include "modules/mapred/path_s3_stub.h"

// this method is defined in .cpp because std::unique_ptr requires a complete type for m_impl
path::path() {}

// this method is defined in .cpp because std::unique_ptr requires a complete type for m_impl
path::~path() {}

// this method is defined in .cpp because std::unique_ptr requires a complete type for m_impl
path::path(const path& rhs)
	: path(rhs.url())
{}

// this method is defined in .cpp because std::unique_ptr requires a complete type for m_impl
path& path::operator=(const path& rhs)
{
	if (rhs.m_impl) {
		m_impl = std::unique_ptr<path_impl>(rhs.m_impl->clone());
	}
	return *this;
}

path::path(std::unique_ptr<path_impl> impl)
	: m_impl(std::move(impl))
{}

path::path(const std::string& url)
{
	if (url == "") {
		return;
	}
	else if (boost::starts_with(url, "s3://")) {
		// SPLOG_P(LOG_DEBUG, "Creating s3 path for %s", url.c_str());
		if (new_path_s3_impl) {
			m_impl = std::unique_ptr<path_impl>(new_path_s3_impl(url));
		}
	}
	else {
		// SPLOG_P(LOG_DEBUG, "Creating standard path for %s", url.c_str());
		m_impl = std::unique_ptr<path_impl>(new path_file_impl(url));
	}
}

bool path::operator<(const path& rhs) const
{
	return url() < rhs.url();
}

bool path::operator==(const path& rhs) const
{
	return url() == rhs.url();
}

std::string path::url() const
{
	if (!valid()) {
		return "";
	}
	return m_impl->url();
}

path::path_type path::type() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->type();
}

std::string path::bare_path() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->bare_path();
}

std::string path_impl::filename() const
{
	std::string::size_type pos = bare_path().rfind('/');
	if (pos == std::string::npos) {
		pos = 0;
	}
	else {
		pos = pos + 1;
	}
	return bare_path().substr(pos);
}

std::unique_ptr<path_impl> path_impl::append(const std::string& suffix) const
{
	std::unique_ptr<path_impl> result(clone());
	auto raw_path = bare_path();
	if (!raw_path.empty() && raw_path[raw_path.size() - 1] != '/') {
		result->m_path += '/';
	}
	result->m_path += suffix;
	return result;
}

std::unique_ptr<path_impl> path_impl::append_unique(const std::string& prefix) const
{
	std::unique_ptr<path_impl> ret;
	do {
		std::string subpath = prefix + "_";
		for (int i = 0; i < 6; i++) {
			subpath += 'a' + (rand() % 26);
		}
		ret = append(subpath);
	} while (ret->exists() != path::e_no_exist);
	return ret;
}

std::string path::filename() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->filename();
}

path path::append(const std::string& suffix) const
{
	if (!valid()) {
		throw io_exception("Appending to empty path");
	}
	return path(m_impl->append(suffix));
}

path path::append_unique(const std::string& prefix) const
{
	if (!valid()) {
		throw io_exception("Appending to empty path");
	}
	return path(m_impl->append_unique(prefix));
}

std::unique_ptr<readable> path::read() const
{
	if (!valid()) {
		throw io_exception("Reading from empty path");
	}
	return m_impl->read();
}

std::unique_ptr<writable> path::write(const path_write_options& options) const
{
	if (!valid()) {
		throw io_exception("Writing to empty path");
	}
	return m_impl->write(options);
}

std::unique_ptr<waiter> path::write_inverted(reset_readable& source, size_t size, const path_write_options& options) const
{
	if (!valid()) {
		throw io_exception("Inverse writing to empty path");
	}
	return m_impl->write_inverted(source, size, options);
}

void path::move(const path& src, const path& dest)
{
	if (!src.valid()) {
		throw io_exception("invalid src path");
	}
	if (!dest.valid()) {
		throw io_exception("invalid dest path");
	}
	if (src.type() != dest.type()) {
		throw io_exception("Cross filesystem type move is not supported");
	}
	src.m_impl->move(src, dest);
}

void path::copy(const path& src, const path& dest, const path_write_options& options)
{
	if (!src.valid()) {
		throw io_exception("invalid src path");
	}
	if (!dest.valid()) {
		throw io_exception("invalid dest path");
	}
	if (src.type() != dest.type()) {
		throw io_exception("Cross filesystem type copy is not supported");
	}
	src.m_impl->copy(src, dest, options);
}

path::exist_enum path::exists() const
{
	if (!valid()) {
		throw io_exception("invalid path object");
	}
	// SPLOG("Checking for existence of '%s'", url().c_str());
	return m_impl->exists();
}


std::time_t path::modify_time() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->modify_time();
}

size_t path::size() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->size();
}

std::vector<std::string> path::list() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->list();
}

void path::mkdir() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	exist_enum e = exists();
	if (e == e_directory) {
		return;
	}
	if (e == e_file) {
		throw io_exception("File of the same name already exists!");
	}
	m_impl->mkdir();
}

bool path::remove() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->rm();
}

bool path::rmdir(bool recursive) const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}

	if (!recursive) {
		return m_impl->rmdir();
	}

	if (!excluded()) {
		walk([](const walk_params& params)
		{
			if (params.state == path::w_file) {
				params.node.remove();
			} else if (params.state == path::w_dir_leave) {
				params.node.rmdir();
			}
		});
		return true;
	}
	SPLOG("path::rmdir> Excluded path: %s", url().c_str());
	return false;
}

bool path::excluded() const
{
	std::string this_url = url();

	// allow paths that are children of 'good_parents'
	std::vector<std::string> good_parents;
	good_parents = Config::instance().get("path_allow_children", good_parents);
	good_parents.push_back(path(CONF_S(storage_root)).url());

	for (const auto& parent : good_parents) {
		if (this_url.find(parent) == 0) {
			return false;
		}
	}

	return true;
}

std::string path::get() const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}

	std::unique_ptr<readable> reader(this->read());
	char buf[64*1024];
	std::string str;

	try {
		while (true) {
			size_t r = reader->read(buf, 64*1024);
			if (r > 0) {
				str += std::string(buf, r);
			}
			if (r < 64*1024) {
				break;
			}
		}
	}
	catch (const io_exception& e) {
		throw io_exception(printstring("in path '%s':\n%s",
			url().c_str(), e.message().c_str()
		));
	}
	return str;
}

void path::put(const std::string& value, const path_write_options& options) const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}

	std::unique_ptr<writable> writer(this->write(options));
	writer->write(value.data(), value.size());
	writer->close();
}

std::string path::str(const exist_enum& e)
{
	static const std::map<exist_enum, std::string> e2s {
		{ e_no_exist,  "does not exist" },
		{ e_file,      "is a file" },
		{ e_directory, "is a directory" }
	};
	const auto& it = e2s.find(e);
	if (it == e2s.cend()) {
		throw io_exception("unknown exist_enum value, please update string conversion function");
	}
	return it->second;
}

// breadth-first walk through the children of this
void path::walk(walker_f fn) const
{
	if (!valid()) {
		throw io_exception("invalid path");
	}
	return m_impl->walk(fn);
}
