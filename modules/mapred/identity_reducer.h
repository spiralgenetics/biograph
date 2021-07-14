
#ifndef __identity_reducer_h__

#include "modules/mapred/reducer.h"

class identity_reducer : public reducer
{
public:
	identity_reducer(const std::string& params) {}
	~identity_reducer() {}

	void start(const std::string& key, kv_sink& out) override { }
	void add_value(const std::string& key, const std::string& value, kv_sink& out) override { out.write(key, value); }
	void end(kv_sink& out) override { }
};

#endif

