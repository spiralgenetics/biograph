#pragma once

#include <string>

class file_wrapper
{
public:
	// Empty buffer
	file_wrapper() = default;

	file_wrapper(const file_wrapper& rhs) = delete;
	file_wrapper& operator=(const file_wrapper&) = delete;
	
	file_wrapper(file_wrapper&& move_me) noexcept
	{
		m_path.swap(move_me.m_path);
		std::swap(m_fd, move_me.m_fd);
	}
	file_wrapper& operator=(file_wrapper&& rhs) noexcept
	{
		m_path.swap(rhs.m_path);
		std::swap(m_fd, rhs.m_fd);
		return *this;
	}

	// See "man -s2 open" for all the allowed flags and modes.
	file_wrapper(const std::string& path, int open_flags, mode_t mode = 0);  

	~file_wrapper() = default;

	bool is_open() const { return m_fd != -1; }
	explicit operator bool() const { return is_open(); }

	void open(const std::string& path, int open_flags, mode_t mode = 0); 

	void close() noexcept;

	const std::string& path() const { return m_path; }
	
	int get_fd() const { return m_fd; }

private:
	std::string m_path;
	int m_fd = -1;
};
