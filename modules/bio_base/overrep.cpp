#include "modules/bio_base/overrep.h"
#include "modules/io/log.h"

static inline uint32_t even_bases(uint64_t kmer) {
	uint64_t ebases = kmer & 0xccccccccccccccccULL;
	uint64_t folded = (ebases | (ebases >> 34));
	return uint32_t(folded);
}

static inline uint32_t odd_bases(uint64_t kmer) {
	uint64_t obases = kmer & 0x3333333333333333ULL;
	uint64_t folded = (obases | (obases >> 30));
	return uint32_t(folded);
}

static inline uint32_t hamming_dist(uint64_t k1, uint64_t k2)
{
	uint64_t diff = k1 ^ k2;
	uint64_t any_diff = (diff & 0x5555555555555555ULL) | ((diff & 0xaaaaaaaaaaaaaaaaULL) >> 1);
	uint32_t count = __builtin_popcount(any_diff >> 32) + __builtin_popcount(any_diff & 0xffffffff);
	return count;
}

overrep_map::overrep_map(size_t kmer_size) 
	: m_kmer_size(kmer_size),
      m_overreps(track_alloc("overreps")),
      m_half0(track_alloc("overreps")),
      m_half1(track_alloc("overreps"))
{
}

overrep_map::~overrep_map() = default;

void overrep_map::add_overrep(const overrep_t& k)
{
	uint32_t elem = uint32_t(m_overreps.size());
	m_overreps.push_back(k);
	m_half0.emplace(even_bases(k.first), elem);
	m_half1.emplace(odd_bases(k.first), elem);
}

bool overrep_map::find_near(const kmer_t& k, overrep_t& out) const
{
	out.second = 0;
	try_side(k, out);
	try_side(rev_comp(k, m_kmer_size), out);
	return out.second > 0;
}

template<class ColType>
class eq_range_looper {
public:
	eq_range_looper(const ColType& col, const typename ColType::key_type& k)
		: m_range(col.equal_range(k)) {}
	typename ColType::const_iterator begin() { return m_range.first; }
	typename ColType::const_iterator end() { return m_range.second; }
private:
	std::pair<typename ColType::const_iterator, typename ColType::const_iterator> m_range;
};

template<class ColType>
eq_range_looper<ColType> make_eq_range(const ColType& col, const typename ColType::key_type& k) {
	return eq_range_looper<ColType>(col, k);
}
	
void overrep_map::try_side(const kmer_t& k, overrep_t& out) const
{
	for(const auto& kvp : make_eq_range(m_half0, even_bases(k))) {
		overrep_t ov = m_overreps[kvp.second];
		if (hamming_dist(k, ov.first) == 1) {
			if (ov.second > out.second) {
				out = ov;
			}
		}
	}
	for(const auto& kvp : make_eq_range(m_half1, odd_bases(k))) {
		overrep_t ov = m_overreps[kvp.second];
		if (hamming_dist(k, ov.first) == 1) {
			if (ov.second > out.second) {
				out = ov;
			}
		}
	}
}

