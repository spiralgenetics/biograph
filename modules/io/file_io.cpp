#include "modules/io/file_io.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"

#include <fstream>
#include <sstream>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

file_reader::file_reader()
	: m_user_file(false)
	, m_file(NULL)
{}

file_reader::file_reader(FILE* file)
	: m_user_file(true)
	, m_file(file)
{}

file_reader::file_reader(const std::string& filename)
	: m_user_file(false)
{
	m_file = fopen(filename.c_str(), "rb");
	if (!m_file) {
		throw io_exception(printstring("%s: in %s", strerror(errno), filename.c_str()));
	}
}

file_reader::file_reader(file_reader&& rhs)
	: m_user_file(rhs.m_user_file)
	, m_file(rhs.m_file)
{
	rhs.m_user_file = false;
	rhs.m_file = NULL;
}

file_reader& file_reader::operator=(file_reader&& rhs)
{
	std::swap(m_user_file, rhs.m_user_file);
	std::swap(m_file, rhs.m_file);
	return *this;
}

file_reader::~file_reader()
{
	close();
}

uint64_t file_reader::size() const
{
	if (m_size >= 0) { return m_size; }

	off_t cur = ftello(m_file);
	if (cur < 0) {
		throw io_exception("Unable to find file pos to get size");
	}
	if (fseeko(m_file, 0, SEEK_END) < 0) {
		throw io_exception("Unable to seek file");
	}
	off_t end = ftello(m_file);
	if (end < 0) {
		throw io_exception("Unable to find file pos to get size (2)");
	}
	if (fseeko(m_file, cur, SEEK_SET) < 0) {
		throw io_exception("Unable to seek file");
	}
	m_size = end;
	return end;
}

void file_reader::seek(uint64_t offset)
{
	if (fseeko(m_file, off_t(offset), SEEK_SET) < 0) {
		throw io_exception("Unable to seek file");
	}
}

size_t file_reader::pos() const
{
	return ftell(m_file);
}

void file_reader::close()
{
	if (m_user_file || m_file == NULL) {
		return;
	}
	if (fclose(m_file) != 0) {
		throw io_exception(strerror(errno));
	}
	m_file = NULL;
}

int file_reader::base_read(char* buf, size_t len)
{
	return fread(buf, 1, len, m_file);
}

file_writer::file_writer(FILE* file)
	: m_user_file(true)
	, m_file(file)
	, m_filename("")
{}

file_writer::file_writer(const std::string& filename, bool append)
	: m_user_file(false),
	m_filename(filename)
{
	m_file = fopen(filename.c_str(), append ? "ab" : "wb");
	if (!m_file) {
		throw io_exception(printstring("Unable to open file %s for writing", filename.c_str()));
	}
}

file_writer::file_writer(file_writer&& rhs)
	: m_user_file(rhs.m_user_file)
	, m_file(rhs.m_file)
	, m_filename(rhs.m_filename)
{
	rhs.m_user_file = false;
	rhs.m_file = NULL;
}

file_writer& file_writer::operator=(file_writer&& rhs)
{
	std::swap(m_user_file, rhs.m_user_file);
	std::swap(m_file, rhs.m_file);
	std::swap(m_filename, rhs.m_filename);
	return *this;
}

file_writer::~file_writer()
{
	close();
}

int file_writer::base_write(const char* buf, int len)
{
	int actually_written = fwrite(buf, 1, len, m_file);
	if (actually_written < len) {
		throw io_exception(printstring("%s: in file %s", strerror(errno), m_filename.c_str()));
	}
	return actually_written;
}

int file_writer::base_flush()
{
	if (m_file == NULL) {
		return 0;
	}
	return fflush(m_file);
}

int file_writer::base_close()
{
	int r = 0;
	if (!m_user_file && m_file != NULL) {
		r = fclose(m_file);
		m_file = NULL;
	}
	return r;
}

void file_writer::seek(uint64_t offset)
{
	if (fseeko(m_file, off_t(offset), SEEK_SET) < 0) {
		throw io_exception("Unable to seek file in writer");
	}
}

size_t file_writer::pos() const
{
	return ftell(m_file);
}

off_t fsize(const std::string& filepath)
{
	struct stat st;
	if (stat(filepath.c_str(), &st) == -1) {
		throw io_exception(printstring("File size error: %s, for path: %s", strerror(errno), filepath.c_str()));
	}
	return st.st_size;
}

std::string slurp_file(const std::string& filepath)
{
	std::ifstream ifs(filepath);
	std::stringstream buf;
	buf << ifs.rdbuf();
	return buf.str();
}
