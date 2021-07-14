#pragma once

#include "modules/io/io.h"
#include "modules/io/progress_tracker.h"

#include <zlib.h>
#include <memory>
#include <array>
#include <future>
#include <mutex>

class zip_reader : public read_wrapper
{
public:
	zip_reader(const zip_reader&) = delete;
	zip_reader(zip_reader&&) = delete;
	zip_reader& operator=(const zip_reader&) = delete;
	zip_reader& operator=(zip_reader&&) = delete;

	zip_reader(readable& source, progress_t& update = no_update);
	~zip_reader();

private:
	int base_read(char* buf, size_t len) override;
	bool check_eof();

  void run_read_thread();
  void close_read_thread();
	int read_internal(char* buf, size_t len);

  static constexpr size_t k_read_buf_size = 64*1024;
  static constexpr size_t k_decompress_buf_size = 64*1024;

	readable& m_source;
	progress_tracker m_tracker;

  std::future<void> m_read_thread;

  std::mutex m_mu;
  std::condition_variable m_read_buffer_consumed;
  std::condition_variable m_read_buffer_avail;
	z_stream m_stream;

  // Set by decompress thread when EOF is encountered.
  bool m_eof = false;

  // Set by main thread when object is destroyed, so that the
  // decompress thread exits even if it's not at EOF.
  bool m_closing = false;

  // Decompressed data that has yet to be read by client.
  char* m_out_buffer = nullptr;
  size_t m_out_size = 0;

	std::array<char, k_read_buf_size> m_buf;
};

// compresses via zlib.
//
// !! Important !!
//
// the GZIP stream will not be terminated (Z_STREAM_END) until after 'close()' is called.
// The destructor of zip_writer will call 'close()'.
// But if you wish to inspect 'sink' before ~zip_writer() is called
// you may need to call 'close()' explicitly to flush the output of the underlying zlib stream.
//
// close() is idempotent.
//
class zip_writer : public write_wrapper
{
public:
	zip_writer(const zip_writer&) = delete;
	zip_writer(zip_writer&&) = delete;
	zip_writer& operator=(const zip_reader&) = delete;
	zip_writer& operator=(zip_writer&&) = delete;

	zip_writer(
		writable& sink,
		progress_t& update_fn = no_update,
		int compression_level = Z_DEFAULT_COMPRESSION,
		int compression_strategy = Z_DEFAULT_STRATEGY
	);
	~zip_writer();

private:
	int base_close() override;
	int base_write(const char* buf, int len) override;

	int compress(int flush);

private:
	writable& m_sink;
	progress_tracker m_tracker;
	z_stream m_stream;
	std::array<char, 64*1024> m_buf;
	bool m_closed = false;
};
