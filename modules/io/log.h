#pragma once

#include <syslog.h>
#include <boost/optional/optional_io.hpp>

// Outputs a message to the current spiral log target.
void spiral_log(int priority, const char* format, ...) __attribute__ ((format (printf, 2, 3)));

// Sets the target for spiral logging messages.
void set_spiral_logging_target(
    const std::function<void(int /* priority */, std::string /* message */)>&);

// log_init sets the spiral logging target to the given file descriptor, or to syslog if fd is -1.
// "name" should identify the current process.
void log_init(const char* name, int fd, bool write_debug = false);
void log_change_name(const char* name);

// Dump build stamp information to the log.
void log_build_stamp();

#define SPLOG(format, ...) \
	spiral_log(LOG_INFO, format, ##__VA_ARGS__);

#define SPLOG_P(priority, format, ...) \
	spiral_log(priority, format, ##__VA_ARGS__);
