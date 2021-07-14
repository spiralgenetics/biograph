#pragma once

#include "modules/io/registry.h"
#include "modules/io/keyvalue.h"
#include "modules/io/log.h"
#include "modules/mapred/task.h"

#include <map>

class reducer 
{
public:
	virtual ~reducer() {}
	virtual void set_watchdog(const std::function<void()>& watchdog) {}
	virtual void start(const std::string& key, kv_sink& context) = 0;
	virtual void add_value(const std::string& key, const std::string& value, kv_sink& context) = 0;
	virtual void end(kv_sink& context) = 0;
	virtual void finalize(kv_sink& context) { }
	virtual std::string get_meta() const { return ""; }
	virtual bool is_summary() { return false; }
	virtual void summarize(std::string& total, const std::string& add) { throw io_exception("Not implemented"); }
	virtual std::string combine_meta(const std::string& m1, const std::string& m2) const { return ""; } 

	virtual task_requirements get_requirements()
	{
		return task_requirements {
			.profile = "normal",
			.cpu_minutes = 10,
		};
	}
};

template<class Subclass, class InKey, class InValue, class OutKey, class OutValue>
class typed_reducer : public reducer 
{
public:
	void output(const OutKey& key, const OutValue& value)
	{
		m_out_context->write_msgpack(key, value);
	}

	void start(const std::string& key_str, kv_sink& context) override
	{
		InKey key;      
		msgpack_deserialize(key, key_str);
		m_out_context = &context;
		((Subclass*) this)->typed_start(key);
		m_out_context = NULL;
	}
	
	void add_value(const std::string& key_str, const std::string& value_str, kv_sink& context) override
	{
		InKey key;      
		InValue value;
		msgpack_deserialize(key, key_str);
		msgpack_deserialize(value, value_str);
		m_out_context = &context;
		((Subclass*) this)->typed_add_value(key, value);
		m_out_context = NULL;
	}
	
	void end(kv_sink& context) override
	{
		m_out_context = &context;
		((Subclass*) this)->typed_end();
		m_out_context = NULL;
	}

protected:
	typed_reducer() {}
	kv_sink* m_out_context = nullptr;
};

template<class Subclass, class Value>
class simple_reducer : public reducer 
{
public:
	void start(const std::string& key, kv_sink& context) override
	{
		m_key = key;
		m_value = Value();
	}

	void add_value(const std::string& key_str, const std::string& value_str, kv_sink& context) override
	{
		Value value;
		msgpack_deserialize(value, value_str);
		((Subclass*) this)->typed_summarize(m_value, value);
	}

	void end(kv_sink& context) override
	{
		context.write(m_key, msgpack_serialize(m_value));
	}

	bool is_summary() override { return true; }

    void summarize(std::string& total, const std::string& add) override
	{
		Value v1;
		Value v2;
		msgpack_deserialize(v1, total);
		msgpack_deserialize(v2, add);
		((Subclass*) this)->typed_summarize(v1, v2);
		total = msgpack_serialize(v1);
	}

private:
	std::string m_key;
	Value m_value;
};

DECLARE_REGISTRY_1(reducer, std::string const&);
