#pragma once

#include "modules/io/transfer_object.h"
#include "json_spirit.h"
#include "modules/io/log.h"

namespace js=json_spirit;

void throw_bad_type(const js::Value_type& found_type, const js::Value_type& expected_type);

template <class T>
js::mValue json_wrap(T& val);

template <class ST, class T>
struct json_wrap_impl
{
	static js::mValue wrap(T& val)
	{
		return js::mValue(transfer_info<T>::get(val));
	}
};

template <class T>
struct json_wrap_impl<transfer_type_null, T>
{
	static js::mValue wrap(T& val) { return js::mValue(); }
};

class json_serialize_context
{
public:
	void set_version(size_t version)
	{
		if (version != 0) {
			m_obj["_ver"] = (uint64_t) version;
		}
		m_version = version;
	}

	bool is_serialize() { return true; }
	bool human_readable() { return true; }
	bool has_field(const std::string& name, size_t tag) { return true; }
	bool is_null(const std::string& name, size_t tag) { return false; }
	size_t get_version() { return m_version; }

	template <class T>
	void transfer_field(const std::string& name, size_t tag, T& obj)
	{
		m_obj[name] = json_wrap(obj);
	}
	const js::mObject& get_result() { return m_obj; }

private:
	size_t m_version = 0;
	js::mObject m_obj;
};

template <>
struct json_wrap_impl<transfer_type_object, js::mValue>
{
	static js::mValue wrap(js::mValue& val)
	{
		return val;
	}
};

template <class T>
struct json_wrap_impl<transfer_type_object, T>
{
	static js::mValue wrap(T& val)
	{
		json_serialize_context ctx;
		transfer_info<T>::object_transfer(ctx, val);
		return ctx.get_result();
	}
};

template <class T>
struct json_wrap_impl<transfer_type_array, T>
{
	static js::mValue wrap(T& val)
	{
		typedef typename transfer_info<T>::iterator iterator;
		iterator itEnd = transfer_info<T>::end(val);
		js::mArray r;
		for(iterator it = transfer_info<T>::begin(val); it != itEnd; ++it) {
			r.push_back(json_wrap(*it));
		}
		return r;
	}
};

template <class T>
struct json_wrap_impl<transfer_type_tuple, T>
{
	static void apply_tuple(js::mArray& r, const boost::tuples::null_type& nothing)
	{}

	template <class Head, class Tail>
	static void apply_tuple(js::mArray& r, boost::tuples::cons<Head, Tail>& x)
	{
		r.push_back(json_wrap(x.get_head()));
		apply_tuple(r, x.get_tail());
	}

	static js::mValue wrap(T& val)
	{
		js::mArray r;
		typename transfer_info<T>::tuple_type tuple = transfer_info<T>::as_boost_tuple(val);
		apply_tuple(r, tuple);
		return r;
	}
};

template <class T>
struct json_wrap_impl<transfer_type_map_object, T>
{
	static js::mValue wrap(T& val)
	{
		js::mObject r;
		typename T::iterator itEnd = val.end();
		for (typename T::iterator it = val.begin(); it != itEnd; ++it) {
			r[it->first] = json_wrap(it->second);
		}
		return r;
	}
};

template <class T>
js::mValue json_wrap(T& val)
{
	typedef typename transfer_info<T>::type sub_type;
	return json_wrap_impl<sub_type, T>::wrap(val);
}

template <class T>
std::string json_serialize(const T& obj, bool pretty = false)
{
	js::mValue json = json_wrap(const_cast<T&>(obj));
	if (pretty) {
		return js::write_formatted(json);
	}
	else {
		return js::write(json);
	}

}

template <class T>
void json_unwrap(T& val, const js::mValue& json);

template <class ST, class T>
struct json_unwrap_impl
{
	static void unwrap(T& val, const js::mValue& json)
	{
		try {
			transfer_info<T>::put(val, json.get_value<ST>());
		}
		catch (const std::runtime_error& e) {
			throw deserialization_error(printstring("Invalid type: %s", e.what()));
		}
	}
};

class json_deserialize_context
{
public:
	json_deserialize_context(const js::mObject& obj)
		: m_obj(obj)
	{}

	void set_version(size_t version) {}
	bool is_serialize() { return false; }
	bool human_readable() { return true; }

	bool has_field(const std::string& name, size_t tag)
	{
		return m_obj.find(name) != m_obj.end();
	}

	bool is_null(const std::string& name, size_t tag)
	{
		return m_obj.find(name)->second.type() == js::null_type;
	}

	size_t get_version()
	{
		if (!has_field("_ver", 0)) {
			return 0;
		}
		uint64_t version;
		json_unwrap(version, m_obj.find("_ver")->second);
		return version;
	}

	template <class T>
	void transfer_field(const std::string& name, size_t tag, T& obj)
	{
		const js::mObject::const_iterator it = m_obj.find(name);
		if (it == m_obj.end()) {
			throw deserialization_error("Missing field: " + name);
		}
		const js::mValue& val = it->second;
		if (val.is_null()) {
			throw deserialization_error("Null field: " + name);
		}
		json_unwrap(obj, val);
	}

private:
	const js::mObject& m_obj;
};

template <class T>
struct json_unwrap_impl<transfer_type_object, T>
{
	static void unwrap(T& val, const js::mValue& json)
	{
		if (json.type() != js::obj_type) {
			throw_bad_type(json.type(), js::obj_type);
		}
		json_deserialize_context ctx(json.get_obj());
		transfer_info<T>::object_transfer(ctx, val);
	}
};

template <>
struct json_unwrap_impl<transfer_type_object, js::mValue>
{
	static void unwrap(js::mValue& val, const js::mValue& json)
	{
		val = json;
	}
};

template <class T>
struct json_unwrap_impl<transfer_type_array, T>
{
	static void unwrap(T& val, const js::mValue& json)
	{
		typedef typename transfer_info<T>::value_type value_type;

		if (json.type() != js::array_type) {
			throw_bad_type(json.type(), js::array_type);
		}

		val = T();

		const js::mArray& arr = json.get_array();
		size_t arr_size = arr.size();
		for (size_t i = 0; i < arr_size; i++) {
			value_type r;
			json_unwrap(r, arr[i]);
			transfer_info<T>::push_back(val, r);
		}
	}
};

template <class T>
struct json_unwrap_impl<transfer_type_tuple, T>
{
	static void apply_tuple(const js::mArray& r, size_t i, const boost::tuples::null_type& nothing)
	{}

	template <class Head, class Tail>
	static void apply_tuple(const js::mArray& r, size_t i, boost::tuples::cons<Head, Tail>& x)
	{
		json_unwrap(x.get_head(), r[i]);
		apply_tuple(r, i+1, x.get_tail());
	}

	static void unwrap(T& val, const js::mValue& json)
	{
		if (json.type() != js::array_type) {
			throw_bad_type(json.type(), js::array_type);
		}
		const js::mArray& arr = json.get_array();
		typedef typename transfer_info<T>::tuple_type tuple_type;
		if (arr.size() != (size_t) boost::tuples::length<tuple_type>::value) {
			throw deserialization_error("Invalid number of tuple elements");
		}
		tuple_type out_tuple;
		apply_tuple(arr, 0, out_tuple);
		transfer_info<T>::from_boost_tuple(val, out_tuple);
	}
};

template <class T>
struct json_unwrap_impl<transfer_type_map_object, T>
{
	static void unwrap(T& val, const js::mValue& json)
	{
		typedef typename transfer_info<T>::value_type value_type;

		if (json.type() != js::obj_type) {
			throw_bad_type(json.type(), js::obj_type);
		}

		val = T();
		const js::mObject& obj = json.get_obj();

		for (const auto& item : obj) {
			value_type item_value;
			json_unwrap(item_value, item.second);
			val[item.first] = item_value;
		}
	}
};

template <class T>
void json_unwrap(T& val, const js::mValue& json)
{
	typedef typename transfer_info<T>::type sub_type;
	json_unwrap_impl<sub_type, T>::unwrap(val, json);
}

template <class T>
void json_deserialize(T& out, const std::string& str)
{
	js::mValue json;
	if (!js::read(str, json)) {
		try {
			js::read_or_throw(str, json);
		}
		catch (const js::Error_position& ep) {
			throw deserialization_error(printstring("Failed on line %u, column %u because %s (%s)",
				ep.line_, ep.column_, ep.reason_.c_str(), str.c_str())
			);
		}
		catch (const std::string& msg) {
			throw deserialization_error(msg);
		}
	}
	json_unwrap(out, json);
}

template <class T>
T inline_json_deserialize(const std::string& str)
{
	T result;
	json_deserialize(result, str);
	return result;
}
