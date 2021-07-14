#pragma once

#include <map>
#include <string>
#include <gtest/gtest.h>

class path;
class manifest;

//generate a manifest at path 'p' that points to 'numKV' unique random key-value pairs.
// The kv pairs are distributed in chunk of size 'chunkSize'.
// Both keys and values are strings of 'kv_size' characters.
// Each character belongs to [a-z].
// A convenience map 'verify' will be filled with the identical pairs used to fill the chunks.
// Pseudo-random generator can be controlled with the seed value. If not provided, this will be to time(NULL).
void gen_random_kv(
	const path& p,
	const size_t numKV,
	const size_t chunkSize,
	const size_t kv_size,
	std::map<std::string, std::string>& verify,
	manifest& out,
	const std::string& encoding,
	unsigned int seed = 0);

// increments a string that is composed of characters in [a-z].
// for instance inc("bg") -> "bh". Handles overflow: inc("az") -> "ba".
std::string inc(const std::string& in);

std::string make_path(const std::string& suffix);

bool diff(const std::string& file1, const std::string& file2);
bool diff(const path& file1, const std::string& file2);
bool zdiff(const std::string& file1, const std::string& file2);
