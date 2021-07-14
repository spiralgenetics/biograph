#include "base/base.h"
#include "modules/io/readable_PRNG.h"
#include <algorithm> 
#include <cassert> 
//#include "io/log.h" 

readable_PRNG::readable_PRNG(size_t size,
		unsigned char randomness,
		unsigned int seed)
	: m_size(size)
	, m_total_read(0)
	, m_randomness(randomness)
	, m_seed(seed)
{
	CHECK_LE( m_randomness, 8 );

	//
	// DR = 2^(randomness) -1
	//
	if (randomness == 0) {
		m_dynamic_range = 0;
	}
	else {
		m_dynamic_range = 1;
		for (unsigned char i=0; i<randomness; i++) {
			m_dynamic_range *=2;
		}
		m_dynamic_range--;
	}

	//SPLOG("dynamic_range=%u", m_dynamic_range);
	CHECK_LE( m_dynamic_range, 255 );
	reset();
}

size_t readable_PRNG::read(char* buf, size_t len)
{
	size_t remainder = m_size - m_total_read;
	size_t written = (remainder > len)? len : remainder;
		
	std::generate_n(buf, written, [&]() { 
		return std::rand() % (m_dynamic_range+1); 
	});

	m_total_read += written;
	return written;
}

void readable_PRNG::reset()
{
	std::srand(m_seed);
	m_total_read = 0;
}

size_t readable_PRNG::size() const
{
	return m_size;
}
