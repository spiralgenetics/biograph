
#pragma once

#include "modules/io/io.h"

class tunstall
{
public:
	static const size_t buf_size(size_t size) { return ((2*size-1) + 7) / 8; }
	// Build a tunstall table of 'size' elements
	tunstall(double one_prob, size_t size);
	// Load a tunstall table from buffer of size bytes
	tunstall(const uint8_t* buf, size_t size);
	// Size of table
	size_t size() { return m_entries.size(); }
	// Write a tunstall table to a buffer, must be buf_size(size()) large
	void write(uint8_t* buf);
	// Return the i'th entry
	const std::vector<bool>& operator[](size_t i) const { return m_entries[i]; }
	// Encode some data
	void encode(std::vector<uint16_t>& out, const uint8_t* buf, size_t buf_size);
	// Decode some data
	void decode(const std::vector<uint16_t>& in, uint8_t* buf, size_t buf_size);
private:
	class bit_writer;
	class bit_reader;
	struct tnode
	{
		tnode* left = nullptr;
		tnode* right = nullptr;
		size_t index = size_t(-1);
	};
	tnode* m_top;
	std::vector<std::vector<bool>> m_entries;
	void make_entries(tnode* cur, std::vector<bool>& bits);
	void write_rec(bit_writer& out, tnode* cur);
	tnode* read_rec(bit_reader& in);
};

