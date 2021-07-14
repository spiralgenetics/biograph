#include <unordered_map>

#include "modules/io/prefix_sum.h"
#include "modules/io/range_coder.h"

// Tracks the probability of a set of options.
class dyn_prob_codec
{
public:
	// Construct a prob_codec with a universe of N symbols
	// On entry of a new symbol, give it a weight of on_first
	dyn_prob_codec(symbol_t universe, uint32_t on_first = 1);
	// Range encode/decode
	void encode(range_encoder& r, symbol_t s);
	void decode(range_decoder& r, symbol_t& s);
	// Update tables given that an option occured
	// called internally, can also be called externally to 'prewarm'
	void update(symbol_t s);
private:
	// Map from symbols to options
	std::unordered_map<uint32_t, uint32_t> m_sym_to_opt;
	// Reverse mapping
	std::vector<uint32_t> m_opt_to_sym;
	// Accumulated option counts + normalization
	prefix_sum_dist m_dist;
	// How much to preweight new entries
	uint32_t m_on_first;
	// Uniform distribution for adding symbols
	uniform_dist m_uniform;
};

class dyn_markov_codec
{
public:
	// Construct a Markov codec, we always begin in state 0
	dyn_markov_codec(symbol_t universe);

	// Range encode/decode, change state
	void encode(range_encoder& r, symbol_t s);
	void decode(range_decoder& r, symbol_t& s);

	// Called internally or to shift state without encoding
	// ie, to prewarm known states
	void update(symbol_t s);

private:
	// find or emplace
	dyn_prob_codec& at(symbol_t s);

	size_t m_universe;
	symbol_t m_cur_state;
	std::unordered_map<symbol_t, dyn_prob_codec> m_states;
};

