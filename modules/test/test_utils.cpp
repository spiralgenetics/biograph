#include <boost/format.hpp>
#include <cstdio>
#include <iostream>
#include <memory>
#include <signal.h>
#include <stdarg.h>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <time.h>

#include "modules/io/config.h"
#include "modules/mapred/base_chunker.h"
#include "modules/mapred/kv_hold.h"
#include "modules/test/test_utils.h"

void gen_random_kv(
	const path& manifest_path,
	const size_t numKV,
	const size_t chunkSize,
	const size_t kv_size,
	std::map<std::string, std::string>& verify,
	manifest& out,
	const std::string& encoding,
	unsigned int seed)
{
	srand(seed ? seed : time(NULL));
	base_chunker<kv_hold> out_chunker("", manifest_path.append("input"), "chunk", chunkSize, 0, out, encoding);
	std::string key;
	std::string value;

	key.resize(kv_size);
	value.resize(kv_size);
	for (size_t i = 0; i < numKV; i++) {
		do {
			for (size_t i = 0; i < kv_size; i++) {
				key[i] = 'a' + random()%26;
				value[i] = 'a' + random()%26;
			}
		} while (verify.find(key) != verify.end());
		verify[key] = value;
		out_chunker.write(key, value);
	}
	out_chunker.close();
}

std::string inc(const std::string& in)
{
	std::string ret = in;
	size_t pos = ret.size() -1;
	char max = 'a'+25;
	 while (true) {
		if (ret[pos] == max) {
			// overflow
			ret.replace(pos, 1, 1, 'a');
			if (pos == 0) {
				break;
			}
			// apply the carry to the next char to the left
			pos--;
		}
		else {
			char last = ret[pos];
			last++;
			ret.replace(pos, 1, 1, last);
			break;
		}
	}
	return ret;
}

std::string make_path(const std::string& suffix)
{
	if (getenv("TEST_TMPDIR")) {
		return std::string(getenv("TEST_TMPDIR")) + "/" + suffix;
	} else {
		// TODO(nils): Remove this case and the "test_root" config
		// when all tests use bazel.
		return CONF_S(test_root) + "/" + suffix;
	}
}

bool diff(const std::string& file1, const std::string& file2)
{
	auto cmd = boost::str(boost::format("diff %s %s") % file1 % file2);
	int sys_ret = system(cmd.c_str());
	int ret = WEXITSTATUS(sys_ret);
	return (ret == 0);
}

bool diff(const path& file1, const std::string& file2)
{
	return	diff(file1.bare_path(), file2);
}

bool zdiff(const std::string& file1, const std::string& file2)
{
	auto cmd = boost::str(boost::format("zdiff %s %s") % file1 % file2);
	int sys_ret = system(cmd.c_str());
	int ret = WEXITSTATUS(sys_ret);
	return (ret == 0);
}
