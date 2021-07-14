#pragma once

#include "modules/io/keyvalue.h"
#include "modules/io/registry.h"
#include "modules/mapred/task.h"
#include "modules/mapred/metadata.h"

class manifest;
class mapper
{
public:
	virtual ~mapper() = default;

	virtual void map(const std::string& key, const std::string& value, kv_sink& context) = 0;
	virtual void set_watchdog(const std::function<void()>& watchdog) {}
	virtual void setup() {}
	virtual void install_metadata(meta::data& metadata) {}

	virtual task_requirements get_requirements()
	{
		return task_requirements {
			.profile = "normal",
			.cpu_minutes = 10,
		};
	}
};

template<class Subclass, class InKey, class InValue, class OutKey, class OutValue>
class typed_mapper : public mapper
{
public:
	void output(const OutKey& key, const OutValue& value)
	{
		m_out_context->write_msgpack(key, value);
	}

	void map(const std::string& key_str, const std::string& value_str, kv_sink& context) override
	{
		InKey key;	
		InValue value;
		msgpack_deserialize(key, key_str);
		msgpack_deserialize(value, value_str);
		m_out_context = &context;
		((Subclass*) this)->typed_map(key, value);
		m_out_context = NULL;
	}
protected:
	typed_mapper() {}

private:
	kv_sink* m_out_context = nullptr;
};

DECLARE_REGISTRY_1(mapper, std::string const&);
