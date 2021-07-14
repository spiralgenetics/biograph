#include <stdio.h>
#include <assert.h>

#include "tools/build_stamp.h"

// We overload build_scm_revision to include both git revision and AWS
// instance ID; see get_workspace_status.
static void get_build_combined_scm_revision(
	std::string& revision, std::string& instance_id) {
	std::string both = BUILD_SCM_REVISION;
	size_t space = both.find(' ');
	assert(space > 0);
	revision = both.substr(0, space);
	instance_id = both.substr(space+1);
}

std::string get_build_scm_revision() {
	std::string revision, instance_id;
	get_build_combined_scm_revision(revision, instance_id);
	return revision;
}

std::string get_build_scm_status() {
	return BUILD_SCM_STATUS;
}

bool build_is_clean() {
	return get_build_scm_status() == "Clean";
}

std::string get_build_instance_id() {
	std::string revision, instance_id;
	get_build_combined_scm_revision(revision, instance_id);
	return instance_id;
}

// AWS gives poor values for current hostname.  So, include the
// instance ID also.
std::string get_build_host() {
	std::string host = BUILD_HOST;

	std::string instance_id = get_build_instance_id();
	if (instance_id.empty()) {
		return host;
	} else if (host == "unknown") {
		return "instance " + instance_id;
	} else {
		return host + " (" + instance_id + ")";
	}
}

std::string get_build_user() {
	return BUILD_USER;
}

time_t get_build_timestamp() {
	// Build timestamp is supplied in seconds since the epoch.
	return BUILD_TIMESTAMP;
}

bool build_info_available() {
	return get_build_timestamp() > 0;
}

