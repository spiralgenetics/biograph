#include "modules/io/version.h"
#include "tools/version.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

// Semantic Versioning 2.0.0
// http://semver.org/

// MAJOR.MINOR.PATCH-PRE+BUILD

// Parse version from a string
product_version::product_version(const std::string& str)
{
	boost::regex re("(\\d+)\\.(\\d+)\\.(\\d+)(?:-([.\\dA-Za-z]+))?(?:\\+([.\\dA-Za-z-]+))?");
	boost::smatch match;
	if (boost::regex_match(str, match, re)) {
		major = atoi(match[1].str().c_str());
		minor = atoi(match[2].str().c_str());
		patch = atoi(match[3].str().c_str());
		pre = match[4].str();
		build = match[5].str();
	} else {
		throw io_exception("Version string does not match regex");
	}
}

std::string product_version::make_string() const
{
	return std::to_string(major)
		+ "." + std::to_string(minor)
		+ "." + std::to_string(patch)
		+ (pre == "" ? "" : ("-" + pre))
		+ (build == "" ? "" : ("+" + build))
	;
}

bool product_version::can_read(const product_version& prev) const
{
	if (major != prev.major) { return false; }
	if (minor < prev.minor) { return false; }
	if (prev.pre != "") {
		if (minor != prev.minor ||
			patch != prev.patch ||
			pre != prev.pre) {
			return false;
		}
	}
	return true;
}

const product_version spec_current_version(SPEC_VERSION);
const product_version seqset_current_version(SEQSET_VERSION);
const product_version biograph_current_version(BIOGRAPH_VERSION);
const product_version biograph_sdk_current_version(BIOGRAPH_VERSION);
