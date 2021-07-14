#pragma once

#include "modules/io/keyvalue.h"
#include "modules/io/registry.h"
#include "modules/mapred/task.h"
#include "modules/mapred/metadata.h"

class dual_mapper
{
public:
	virtual ~dual_mapper() {}
	virtual void dual_map(const std::string& key, const std::string& value, kv_sink& cxt1, kv_sink& cxt2) = 0;
	virtual void set_watchdog(const std::function<void()>& watchdog) {}
	virtual void setup() {}
	virtual void install_metadata1(meta::data& metadata) {}
	virtual void install_metadata2(meta::data& metadata) {}

	virtual task_requirements get_requirements()
	{
		return task_requirements {
			.profile = "normal",
			.cpu_minutes = 1,
		};
	}
};

template<class Subclass, 
	class InKey, class InValue, 
	class OutKey1, class OutValue1, 
	class OutKey2, class OutValue2>
class typed_dual_mapper : public dual_mapper
{
public:
	void output1(const OutKey1& key, const OutValue1& value)
	{
		m_out_ctx1->write_msgpack(key, value);
	}

	void output2(const OutKey2& key, const OutValue2& value)
	{
		m_out_ctx2->write_msgpack(key, value);
	}

	void dual_map(const std::string& key_str, const std::string& value_str, kv_sink& ctx1, kv_sink& ctx2) override
	{
		InKey key;
		InValue value;
		msgpack_deserialize(key, key_str);
		msgpack_deserialize(value, value_str);
		m_out_ctx1 = &ctx1;
		m_out_ctx2 = &ctx2;
		((Subclass*) this)->typed_map(key, value);
		m_out_ctx1 = nullptr;
		m_out_ctx2 = nullptr;
	}
private:
	kv_sink* m_out_ctx1 = nullptr;
	kv_sink* m_out_ctx2 = nullptr;
};

DECLARE_REGISTRY_1(dual_mapper, std::string const&);
