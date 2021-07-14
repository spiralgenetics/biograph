#ifndef __readable_PRNG_h__
#define __readable_PRNG_h__

#include "modules/io/io.h"


// Generate un-buffered pseudo-random data.
// Because its data is generated on-the-fly, the memory footprint is constant and tiny,
// regardless of the 'size'.
//
// 'size' is the maximum number of bytes to be generated
//
// 'randomness' is in [0..8]
//
// 8 produces data as random as the underlying PRNG, so basically uncompressible data.
// 0 produces constant data equal to 0, so most compressible.
//
// Note that randomness, as defined by this class, is roughly linear:
// randomness(2) is 2 times more compressible than randomness(4).
// This is achieved by constraining the dynamic range of the output of the PRNG.
//
// 'seed' is used to initialize the PRNG
//
// 'reset()' can be called to re-generate the same data
//
// Implementation details:
// std::rand is used for the PRNG. This is not super-random, but is enough to fool zlib.
//
class readable_PRNG : public reset_readable
{
	public:
		readable_PRNG(size_t size, unsigned char randomness, unsigned int seed);
		size_t read(char* buf, size_t len) override;
		void reset() override;
		size_t size() const;

	private:
		size_t m_size;
		size_t m_total_read;
		unsigned char m_randomness;
		unsigned char m_dynamic_range;
		unsigned int m_seed;
};

#endif
