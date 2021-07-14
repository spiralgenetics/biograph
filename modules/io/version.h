#pragma once

#include <string>
#include "modules/io/transfer_object.h"

struct product_version
{
	unsigned int major;
	unsigned int minor;
	unsigned int patch;
	std::string pre;
	std::string build;

	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(major);
		FIELD(minor);
		FIELD(patch);
		FIELD(pre);
		FIELD(build);
	}

	// Empty for deserialization
	product_version() = default;
	// Parse version from a string
	product_version(const std::string& str);
	// Convert version back to a string
	std::string make_string() const;
	// Determine if a given build can read a file from another build
	bool can_read(const product_version& prev) const;
};

inline bool operator==(const product_version& lhs, const product_version& rhs)
{
	return lhs.major == rhs.major
		&& lhs.minor == rhs.minor
		&& lhs.patch == rhs.patch
	;
}

extern const product_version spec_current_version;
extern const product_version seqset_current_version;
extern const product_version biograph_current_version;
extern const product_version biograph_sdk_current_version;
