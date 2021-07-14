
#include "modules/io/dynamic_codecs.h"

dyn_prob_codec::dyn_prob_codec(symbol_t universe, uint32_t on_first)
	: m_dist()
	, m_on_first(on_first)
	, m_uniform(universe)
{
	m_dist.inner().push_back(1);
}

void dyn_prob_codec::encode(range_encoder& r, symbol_t symbol)
{
	// SPLOG("Encoding %d\n", (int) symbol);
	// Search for symbol
	auto it = m_sym_to_opt.find(symbol);
	if (it == m_sym_to_opt.end()) {
		// If it's not already on my list
		if (! m_opt_to_sym.empty()) {
			// SPLOG("Sending out an 'new_opt'\n");
			r.encode(m_dist, 0);
		}
		// Encode 'uniform' version
		// SPLOG("Sending out new symbol\n");
		r.encode(m_uniform, symbol);
	} else {
		// It's on the list, just encode as option #
		// SPLOG("Encoding as option %d\n", (int) it->second);
		r.encode(m_dist, it->second + 1);
	}
	// Add to table
	update(symbol);
}

void dyn_prob_codec::decode(range_decoder& r, symbol_t& symbol)
{
	// Set default to 'new_opt'
	symbol_t option = 0;
	// If not empty, pull option
	if (! m_opt_to_sym.empty()) {
		// decode option # (or new_opt)
		// SPLOG("Decoding option -> %d", (int) option);
		option = r.decode(m_dist);
	}
	if (option == 0) {
		// Was an new_opt, grab uniform
		symbol = r.decode(m_uniform);
		// SPLOG("Got new, decoding uniform -> %d", (int) symbol);
	} else {
		// Get actual symbol
		symbol = m_opt_to_sym[option - 1];
		// SPLOG("Got option, converting to symbol -> %d", (int) symbol);
	}
	// Add to table
	update(symbol);
}

void dyn_prob_codec::update(symbol_t symbol)
{
	// SPLOG("Doing update for symbol %d, new accum:", (int) symbol);
	auto it = m_sym_to_opt.find(symbol);
	if (it == m_sym_to_opt.end()) {
		// Add a new entry
		it = m_sym_to_opt.insert(std::make_pair(symbol, m_opt_to_sym.size())).first;
		m_opt_to_sym.push_back(symbol);
		m_dist.inner().push_back(m_on_first);
	} else {
		m_dist.inner().add(it->second + 1, 1);
	}
}

dyn_markov_codec::dyn_markov_codec(symbol_t universe)
	: m_universe(universe)
	, m_cur_state(0)
{
	m_states.emplace(m_cur_state, m_universe);
}

void dyn_markov_codec::encode(range_encoder& r, symbol_t s)
{
	at(m_cur_state).encode(r, s);
	m_cur_state = s;
}

void dyn_markov_codec::decode(range_decoder& r, symbol_t& s)
{
	at(m_cur_state).decode(r, s);
	m_cur_state = s;
}

void dyn_markov_codec::update(symbol_t s)
{
	at(m_cur_state).update(s);
	m_cur_state = s;
}

dyn_prob_codec& dyn_markov_codec::at(symbol_t s)
{
	return m_states.emplace(s, m_universe).first->second;
}



