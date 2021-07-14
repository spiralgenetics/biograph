#pragma once

#include "modules/io/json_transfer.h"
#include "modules/io/simple_metadata.h"

#include <functional> 
#include <map>
#include <string> 

namespace meta {

namespace ns {

// spiral's internal stuff. Invisible to end-user.
extern const char* internal;

// visible, but read-only to end-user
extern const char* readonly;

// ns for the end-user, R/W.
extern const char* user;

} // namespace ns

namespace merge {

struct init {};
struct params
{
	std::string ns;
	std::string key;
	js::mValue value1;
	js::mValue value2;
};

// prototype for a metadata merge handler.  
// Called when the namespace and keys match for two manifests.
// Returns a std::string which is the resolved value that should be used for (ns, key)
typedef std::function<js::mValue (const params& params)> function;

init register_fn(const std::string& key, function fn);

js::mValue first(const params& params);
js::mValue second(const params& params);
js::mValue sum(const params& params);
js::mValue collide(const params& params);

} // namespace merge

class data : public simple_metadata
{
public:
	typedef std::map<std::string, js::mValue>  key_value_t;
	typedef std::map<std::string, key_value_t> namespaces_t;

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(m_data); 
	}

	bool has_key(const std::string& ns, const std::string& key) const
	{
		return m_data.count(ns) && m_data.at(ns).count(key);
	}

	template <class V> 
	V get(const std::string& ns, const std::string& key) const
	{
		return get_(ns, key, V(), true);
	}

	template <class V>
	V get(const std::string& ns, const std::string& key, const V& default_value) const
	{
		return get_(ns, key, default_value, false);
	}

	template <class V>
	void set(const std::string& ns, const std::string& key, const V& value)
	{
		m_data[ns][key] = json_wrap(const_cast<V&>(value));
	}

	void set(const std::string& ns, const std::string& key, const char* value)
	{
		std::string str(value);
		m_data[ns][key] = json_wrap(str);
	}

  // For simple_metadata where we don't want to pull in all of
  // //modules/mapred if we don't have to.
  void set_simple_json(const std::string& key, const js::mValue& value) override {
    m_data[meta::ns::readonly][key] = value;
  }

	void unset(const std::string& ns, const std::string& key);

	void merge(const data& in);

	const meta::data::namespaces_t& raw() const
	{
		return m_data; 
	}

	// add ('spiral', 'created', <time in ISO 8601 format>) and ('internal', 'created', std::time_t)
	void set_creation_time_now();

	// add command line parameters JSON string to metadata
	void set_options(const std::string& step_name, const std::string& the_options_json);

	// add step runtime to metadata.
	void set_runtime(time_t start_time);

private:
	template <class V> 
	V get_(
		const std::string& ns,
		const std::string& key,
		const V&           default_value,
		const bool         throw_on_error) const
	{
		if (ns.empty()) {
			throw io_exception("empty metadata namespace: <empty>/"+key);
		}

		if (key.empty()) {
			throw io_exception("empty metadata key: "+ns+"/<empty>");
		}

		auto ns_it = m_data.find(ns);
		if (ns_it == m_data.end()) {
			if (throw_on_error) {
				throw io_exception("metadata namespace "+ns+"/"+key+" does not exist");
			}
			return default_value;
		}

		auto key_it = ns_it->second.find(key);
		if (key_it == ns_it->second.cend()) {
			if (throw_on_error) {
				throw io_exception("metadata key "+ns+"/"+key+" does not exist");
			}
			return default_value;
		}

		V result;
		json_unwrap(result, key_it->second);
		return result;
	}

private:
	namespaces_t m_data;
};

} // namespace meta
