#pragma once

#include "modules/mapred/path.h"

struct path_impl
{
	virtual ~path_impl() {}

	virtual path::path_type type() const = 0;
	virtual path_impl* clone() const = 0;
	virtual std::string url() const = 0;
	virtual std::string bare_path() const { return m_path; }
	virtual std::string filename() const;
	virtual std::unique_ptr<path_impl> append(const std::string& suffix) const;
	virtual std::unique_ptr<path_impl> append_unique(const std::string& prefix) const;
	virtual std::unique_ptr<readable> read() const = 0;
	virtual std::unique_ptr<writable> write(const path_write_options& options) const = 0;
	virtual std::unique_ptr<waiter> write_inverted(reset_readable& source, size_t size, const path_write_options& options) const = 0;
	virtual void move(const path& src, const path& dest) const = 0;
	virtual void copy(const path& src, const path& dest, const path_write_options& options) const = 0;
	virtual path::exist_enum exists() const = 0;
	virtual std::time_t modify_time() const = 0;
	virtual size_t size() const = 0;
	virtual std::vector<std::string> list() const = 0;
	virtual void mkdir() const = 0;
	virtual bool rm() const = 0;
	virtual bool rmdir() const = 0;
	virtual void walk(path::walker_f fn) const = 0;

	// this is needed for subclasses (concrete impls) that might want access to other impls
	// this works because path_impl is a friend of path, however not the subclasses
	static path_impl* get_impl(const path& path) { return path.m_impl.get(); }

	std::string m_path;
};

struct path_file_impl : public path_impl
{
	path_file_impl(const std::string& url);

	path::path_type type() const override { return path::FILE; }
	path_impl* clone() const override;
	std::string url() const override;
	std::unique_ptr<readable> read() const override;
	std::unique_ptr<writable> write(const path_write_options& options) const override;
	std::unique_ptr<waiter> write_inverted(reset_readable& source, size_t size, const path_write_options& options) const override;
	void move(const path& src, const path& dest) const override;
	void copy(const path& src, const path& dest, const path_write_options& options) const override;
	path::exist_enum exists() const override;
	std::time_t modify_time() const override;
	size_t size() const override;
	std::vector<std::string> list() const override;
	void mkdir() const override;
	bool rm() const override;
	bool rmdir() const override;
	void walk(path::walker_f fn) const override;
};

