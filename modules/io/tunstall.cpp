
#include "base/base.h"
#include "modules/io/tunstall.h"
#include "modules/io/log.h"
#include <map>
#include <math.h>

class tunstall::bit_writer
{
public:
	bit_writer(uint8_t* bits) : m_next(bits) {}
	~bit_writer() {
		if (m_count == 0) { return; }
		while(m_count != 0) { write(false); }
	}
	void write(bool b) {
		m_cur <<= 1;
		m_cur |= b;
		m_count++;
		if (m_count == 8) {
			*m_next = m_cur;
			m_next++;
			m_cur = 0;
			m_count = 0;
		}
	}
			
private:
	uint8_t* m_next;
	uint8_t m_count = 0;
	uint8_t m_cur = 0;
};

class tunstall::bit_reader
{
public:
	bit_reader(const uint8_t* bits, size_t max_size) : m_next(bits), m_size(max_size) {}
	bool read() {
		if (m_count == 0) {
			if (m_size == 0) {
				throw io_exception("Read off the end of a bit_reader");
			}
			m_cur = *m_next;
			m_next++;
			m_size--;
			m_count = 8;
		}
		m_count--;
		return m_cur & (1 << m_count);
	}
private:
	const uint8_t* m_next;
	size_t m_size;
	uint8_t m_cur = 0;
	uint8_t m_count = 0;
};

tunstall::tunstall(double one_prob, size_t size)
{
	CHECK_GE(size, 2); 

        double zero_prob = 1 - one_prob;
        double one_ent = -log2(one_prob);
        double zero_ent = -log2(zero_prob);

	std::multimap<double, tnode*> leaves;
	m_top = new tnode();
	leaves.emplace(0.0, m_top);
	while(leaves.size() < size) {
		double ent = leaves.begin()->first;
		tnode* n = leaves.begin()->second;
		leaves.erase(leaves.begin());
		n->left = new tnode;
		n->right = new tnode;
		leaves.emplace(ent + zero_ent, n->left);
		leaves.emplace(ent + one_ent, n->right);
	}
	std::vector<bool> empty;
	make_entries(m_top, empty);
}

tunstall::tunstall(const uint8_t* buf, size_t size)
{
	bit_reader br(buf, size);
	m_top = read_rec(br);
	std::vector<bool> empty;
	make_entries(m_top, empty);
}

void tunstall::write(uint8_t* buf)
{
	bit_writer bw(buf);
	write_rec(bw, m_top);
}

void tunstall::encode(std::vector<uint16_t>& out, const uint8_t* buf, size_t buf_size)
{
	bit_reader br(buf, buf_size);
	tnode* cur = m_top;
	for(size_t i = 0; i < buf_size*8; i++) {
		if (br.read()) {
			cur = cur->right;
		} else {
			cur = cur->left;
		}
		if (!cur->left) {
			out.push_back(cur->index);
			cur = m_top;
		}
	}
	if (cur == m_top) { return; }
	while(cur->left) { cur = cur->left; }
	out.push_back(cur->index);
}

void tunstall::decode(const std::vector<uint16_t>& in, uint8_t* buf, size_t buf_size)
{
	bit_writer bw(buf);
	buf_size *= 8;
	for(size_t i = 0; i < in.size(); i++) {
		const std::vector<bool> bits = m_entries[in[i]];
		for(size_t j = 0; j < bits.size(); j++) {
			if (buf_size == 0) {
				return;
			}
			bw.write(bits[j]);
			buf_size--;
		}
	}	
}

void tunstall::make_entries(tnode* cur, std::vector<bool>& bits)
{
	if (!cur->left) {
		cur->index = m_entries.size();
		m_entries.push_back(bits);
		return;
	}
	bits.push_back(0);
	make_entries(cur->left, bits);
	bits.pop_back();
	bits.push_back(1);
	make_entries(cur->right, bits);
	bits.pop_back();
}

void tunstall::write_rec(bit_writer& out, tnode* cur)
{
	if (!cur->left) {
		out.write(false);
		return;
	}
	out.write(true);
	write_rec(out, cur->left);
	write_rec(out, cur->right);
}

tunstall::tnode* tunstall::read_rec(bit_reader& in)
{
	tnode* r = new tnode;
	if (in.read()) {
		r->left = read_rec(in);
		r->right = read_rec(in);
	}
	return r;
}

