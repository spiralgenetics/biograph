#include "modules/io/log.h"
#include <stddef.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "tools/build_stamp.h"

static char log_name[256] = {0};
static int log_fd = -1;
static int log_level = LOG_DEBUG;

static std::function<void(int, const std::string&)> g_log_output_function;

void set_spiral_logging_target(
    const std::function<void(int /* priority */, std::string /* message */)>& new_target) {
  g_log_output_function = new_target;
}

void spiral_log(int priority, const char* format, ...) {
	char log_buffer[PIPE_BUF] = {0};
	va_list args;
	va_start(args, format);
    if (g_log_output_function) {
      char* buf;
      int len = vasprintf(&buf, format, args);
      std::string stringified(buf, len);
      free(buf);
      g_log_output_function(priority, stringified);
    } else if (log_fd >= 0 && priority <= log_level) {
		time_t t = time(NULL);
		char *time_str = ctime(&t);
		// The following line is terrible code!
		time_str[strlen(time_str) - 1] = 0;
		int header = snprintf(log_buffer, PIPE_BUF, "%s %s[%d]: ", time_str, log_name, getpid());
		int rest = vsnprintf(log_buffer + header, PIPE_BUF - header, format, args);
		int total = header + rest;
		if (total >= PIPE_BUF) {
			total = PIPE_BUF - 1;
		}
		log_buffer[total++] = '\n';
		// Avoid pedantic warning
		if (write(log_fd, log_buffer, total) < 0) {}
	} else {
		vsyslog(priority, format, args);
	}
	va_end(args);
}

void log_init(const char* name, int fd, bool write_debug)
{
	log_fd = fd;

	// If debug logging is disabled (default), only log up to LOG_INFO
	if (! write_debug) {
		setlogmask(LOG_UPTO(LOG_INFO));
		log_level = LOG_INFO;
	}

	if (getenv("SPIRAL_LOG_STDERR")) {
		// Redirect all logging to standard error for
		// testing/debugging.
		log_fd = 2;
	}

	if (log_fd == -1) {
		// If no log file specified, log to syslog.
		int option = LOG_PID;
		openlog(name, option, LOG_LOCAL0);
	} else {
		if (name) {
			strncpy(log_name, name, 256);
			log_name[255] = 0;
		}
	}
}

void log_change_name(const char* name)
{
	strncpy(log_name, name, 256);
	log_name[255] = 0;
}

void log_build_stamp()
{
	if (build_info_available()) {
		time_t built_at = get_build_timestamp();
		char *time_str = ctime(&built_at);
		SPLOG("Built at %s by %s on %s from revision %s (%s build)",
			time_str, get_build_user().c_str(), get_build_host().c_str(),
			get_build_scm_revision().c_str(), get_build_scm_status().c_str());
	} else {
		SPLOG("Unversioned development build; DO NOT RELEASE.");
	}
}
