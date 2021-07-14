#pragma once

#include <string>
#include <ctime>

// True if this build contains a build stamp.
bool build_info_available();

// Git commit ID.
std::string get_build_scm_revision();

// False if any files have been modified from the git commit returned
// by get_build_scm_revision().
bool build_is_clean();

// "Clean" if build_is_clean(); "Modified" otherwise.
std::string get_build_scm_status();

// The host used to make this build, or "unknown".
std::string get_build_host();

// The AWS instance ID used to make this build, if built on AWS.
std::string get_build_instance_id();

// The user that executed this build.
std::string get_build_user();

// The time this build was run.
time_t get_build_timestamp();
