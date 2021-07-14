
#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/mapred/reducer.h"

class prefix_reducer : public reducer
{
public:
	prefix_reducer(const std::string& params) {}

	void start(const std::string& key_str, kv_sink& context) override {}
        void add_value(const std::string& key_str, const std::string& value_str, kv_sink& context) override { 
		m_last_key = key_str;
		m_last_value = value_str;
	}
        void end(kv_sink& context) override {
		context.write(m_last_key, m_last_value);
	}

private:
	std::string m_last_key;
	std::string m_last_value;
};

