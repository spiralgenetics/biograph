
#ifndef __msgpack_transfer_h__
#define __msgpack_transfer_h__

#include "modules/io/utils.h"
#include "modules/io/transfer_object.h"

// temporary define for _MSC_VER
// we do this because msgpack needs _MSC_VER to compile
// but defining _MSC_VER conflicts with BOOST+MinGW
#ifdef _WIN32
	#define _MSC_VER 1600
		#include <msgpack.h>
	#undef _MSC_VER
#else
	#include <msgpack.h>
#endif

class msgpack_serialize_context;

template<class T>
void msgpack_wrap(msgpack_packer* pk, T& val);

template<class ST, class T>
struct msgpack_wrap_impl;

#define BASE_TYPE_PACK(transfer_type, packfunc) \
template<class T> \
struct msgpack_wrap_impl<transfer_type, T> \
{ \
	static void wrap(msgpack_packer* pk, T& val) \
	{ \
		packfunc(pk, (transfer_type) transfer_info<T>::get(val)); \
	} \
}; 

inline void msgpack_pack_bool(msgpack_packer* pk, bool val)
{
	if (val)
		msgpack_pack_true(pk);
	else
		msgpack_pack_false(pk);
}

BASE_TYPE_PACK(int64_t, msgpack_pack_int64)
BASE_TYPE_PACK(uint64_t, msgpack_pack_uint64)
BASE_TYPE_PACK(double, msgpack_pack_double)
BASE_TYPE_PACK(bool, msgpack_pack_bool)

template<class T>
struct msgpack_wrap_impl<transfer_type_null, T>
{
	static void wrap(msgpack_packer*pk, T& val) { msgpack_pack_nil(pk); }
};

template<class T>
struct msgpack_wrap_impl<transfer_type_string, T>
{
	static void wrap(msgpack_packer* pk, T& val)
	{
		std::string str = transfer_info<T>::get(val);
		msgpack_pack_raw(pk, str.size());
		msgpack_pack_raw_body(pk, str.data(), str.size());
	}
}; 

struct precheck_context 
{
	precheck_context() : m_count(0), m_version(0) {}
	void set_version(size_t version) { m_version = version; }
	bool is_serialize() { return true; }
	bool human_readable() { return true; }
	bool has_field(const std::string& name, size_t tag) { return true; }
	bool is_null(const std::string& name, size_t tag) { return false; }
	size_t get_version() { return m_version; }
	template<class T> void transfer_field(const std::string& name, size_t tag, T& obj) 
	{ m_count++; }
	
	int m_count;
	size_t m_version;
};

class msgpack_serialize_context
{
public:
	msgpack_serialize_context(msgpack_packer* pk, const precheck_context& pcc) 
		: m_pk(pk)
		, m_version(pcc.m_version) 
	{
		if (m_version != 0)
		{
			msgpack_pack_map(m_pk, pcc.m_count + 1);
			msgpack_pack_nil(m_pk);
			msgpack_pack_uint64(m_pk, pcc.m_version);	
		}
		else
			msgpack_pack_map(m_pk, pcc.m_count);
	
	}
	void set_version(size_t version) {}
	bool is_serialize() { return true; }
	bool human_readable() { return false; }
	bool has_field(const std::string& name, size_t tag) { return true; }
	bool is_null(const std::string& name, size_t tag) { return false; }
	size_t get_version() { return m_version; }
	template<class T>
	void transfer_field(const std::string& name, size_t tag, T& obj)
	{
		msgpack_pack_uint64(m_pk, tag);
		msgpack_wrap(m_pk, obj);
	}

private:
	msgpack_packer* m_pk;
	size_t m_version;
};

template<class T>
struct msgpack_wrap_impl<transfer_type_object, T>
{
	static void wrap(msgpack_packer* pk, T& val)
	{
		precheck_context pcc;
		transfer_info<T>::object_transfer(pcc, val);
		msgpack_serialize_context ctx(pk, pcc);
		transfer_info<T>::object_transfer(ctx, val);
	}
};

template<class T>
struct msgpack_wrap_impl<transfer_type_array, T>
{
	static void wrap(msgpack_packer* pk, T& val)
	{
		size_t count = transfer_info<T>::size(val);
		typedef typename transfer_info<T>::iterator iterator;
		iterator itEnd = transfer_info<T>::end(val);
		msgpack_pack_array(pk, count);
		for(iterator it = transfer_info<T>::begin(val); it != itEnd; ++it)
			msgpack_wrap(pk, *it);  
	}
};

template<class T>
struct msgpack_wrap_impl<transfer_type_tuple, T>
{
	static void apply_tuple(msgpack_packer* pk, const boost::tuples::null_type& nothing) 
	{}

	template <class Head, class Tail>
	static void apply_tuple(msgpack_packer* pk, boost::tuples::cons<Head, Tail>& x) 
	{
		msgpack_wrap(pk, x.get_head());
		apply_tuple(pk, x.get_tail());
	}

	static void wrap(msgpack_packer* pk, T& val)
	{
		typedef typename transfer_info<T>::tuple_type tt;
		size_t length = boost::tuples::length<tt>::value;
		tt tuple = transfer_info<T>::as_boost_tuple(val);
		msgpack_pack_array(pk, length);
		apply_tuple(pk, tuple);
	}
};

template<class T>
struct msgpack_wrap_impl<transfer_type_map_object, T>
{
	static void wrap(msgpack_packer* pk, T& val)
	{
		msgpack_pack_map(pk, val.size());
		typename T::iterator itEnd = val.end();
		for(typename T::iterator it = val.begin(); it != itEnd; ++it)
		{
			msgpack_wrap(pk, const_cast<std::string&>(it->first));
			msgpack_wrap(pk, it->second);
		}
	}
};

template<class T>
void msgpack_wrap(msgpack_packer* pk, T& val)
{
	typedef typename transfer_info<T>::type sub_type;
	msgpack_wrap_impl<sub_type, T>::wrap(pk, val);
}

template<class T>
std::string msgpack_serialize(const T& obj)
{
	msgpack_sbuffer* buffer = msgpack_sbuffer_new();
	msgpack_packer* pk = msgpack_packer_new(buffer, msgpack_sbuffer_write);
	msgpack_wrap(pk, const_cast<T&>(obj));
	std::string r(buffer->data, buffer->size);
	msgpack_packer_free(pk);
	msgpack_sbuffer_free(buffer);
	return r;
}

template<class T>
std::string msgpack_serialize(int type, const T& obj)
{
	msgpack_sbuffer* buffer = msgpack_sbuffer_new();
	msgpack_packer* pk = msgpack_packer_new(buffer, msgpack_sbuffer_write);
	msgpack_pack_array(pk, 2);
	msgpack_pack_int64(pk, type);	
	msgpack_wrap(pk, const_cast<T&>(obj));
	std::string r(buffer->data, buffer->size);
	msgpack_packer_free(pk);
	msgpack_sbuffer_free(buffer);
	return r;
}

inline void printpack(const std::string& buf)
{
	msgpack_unpacked msg;
	msgpack_unpacked_init(&msg);
	msgpack_unpack_next(&msg, buf.data(), buf.size(), NULL);
	msgpack_object_print(stdout, msg.data);
	msgpack_unpacked_destroy(&msg);
}

template<class T>
void msgpack_unwrap(T& val, const msgpack_object& mpo);

template<class ST, class T>
struct msgpack_unwrap_impl;

#define BASE_TYPE_UNPACK(transfer_type, unpacktype, unpackfield) \
template<class T> \
struct msgpack_unwrap_impl<transfer_type, T> \
{ \
	static void unwrap(T& val, const msgpack_object& mpo) \
	{ \
		if (mpo.type == unpacktype) \
			transfer_info<T>::put(val, mpo.via.unpackfield); \
		else \
			throw deserialization_error("msgpack object is not of type "#unpacktype); \
	} \
}; 

BASE_TYPE_UNPACK(transfer_type_boolean, MSGPACK_OBJECT_BOOLEAN, boolean);
BASE_TYPE_UNPACK(transfer_type_real, MSGPACK_OBJECT_DOUBLE, dec);
BASE_TYPE_UNPACK(transfer_type_unsigned, MSGPACK_OBJECT_POSITIVE_INTEGER, u64);
// This one is a bit of an evil misuse of macros, but saves some typing
BASE_TYPE_UNPACK(transfer_type_signed, MSGPACK_OBJECT_POSITIVE_INTEGER || mpo.type == MSGPACK_OBJECT_NEGATIVE_INTEGER, i64);

template<class T>
struct msgpack_unwrap_impl<transfer_type_string, T>
{
	static void unwrap(T& val, const msgpack_object& mpo)
	{
		if (mpo.type != MSGPACK_OBJECT_RAW)
			throw deserialization_error("msgpack object is not of type MSGPACK_OBJECT_RAW");
		std::string s(mpo.via.raw.ptr, mpo.via.raw.size); 
		transfer_info<T>::put(val, s);
	}
}; 

class msgpack_deserialize_context
{
public:
	msgpack_deserialize_context(const msgpack_object_map& map) 
		: m_map(map)
		, m_version(0)
	{
		for(size_t i = 0; i < map.size; i++)
		{
			if (m_map.ptr[i].key.type == MSGPACK_OBJECT_NIL)
			{
				if (m_map.ptr[i].val.type != MSGPACK_OBJECT_POSITIVE_INTEGER)
					throw deserialization_error("Non numeric version #");
				m_version = m_map.ptr[i].val.via.u64;
			}
			else if (m_map.ptr[i].key.type == MSGPACK_OBJECT_POSITIVE_INTEGER)
				m_tag_map[m_map.ptr[i].key.via.u64] = i;
			else
				throw deserialization_error("Map has non-numeric tag");
		}
	}
	void set_version(size_t version) {}
	bool is_serialize() { return false; }
	bool human_readable() { return false; }
	bool has_field(const std::string& name, size_t tag) 
	{
		return m_tag_map.find(tag) != m_tag_map.end();
	}
	bool is_null(const std::string& name, size_t tag) 
	{
		size_t i = m_tag_map.find(tag)->second;
		return m_map.ptr[i].val.type == MSGPACK_OBJECT_NIL;
	}
	size_t get_version() { return m_version; }

	template<class T>
	void transfer_field(const std::string& name, size_t tag, T& obj)
	{
		std::map<size_t, size_t>::iterator it = m_tag_map.find(tag);
		if (it == m_tag_map.end())
			throw deserialization_error("Missing field: " + name);
		size_t i= it->second;
		msgpack_unwrap(obj, m_map.ptr[i].val);
	}

private:
	const msgpack_object_map& m_map;
	size_t m_version;
	std::map<size_t, size_t> m_tag_map;  // Tag -> entry #
};

template<class T>
struct msgpack_unwrap_impl<transfer_type_object, T>
{
	static void unwrap(T& val, const msgpack_object& mpo)
	{
		if (mpo.type != MSGPACK_OBJECT_MAP)
			throw deserialization_error("msgpack object is not of type MSGPACK_OBJECT_MAP");
		msgpack_deserialize_context ctx(mpo.via.map);
		transfer_info<T>::object_transfer(ctx, val);
	}
};

template<class T>
struct msgpack_unwrap_impl<transfer_type_array, T>
{
	static void unwrap(T& val, const msgpack_object& mpo)
	{
		typedef typename transfer_info<T>::value_type value_type;

		if (mpo.type != MSGPACK_OBJECT_ARRAY)
			throw deserialization_error("msgpack object is not of type MSGPACK_OBJECT_ARRAY");

		val = T();
		const msgpack_object_array& arr = mpo.via.array;
		for(size_t i = 0; i < arr.size; i++)
		{
			value_type r;
			msgpack_unwrap(r, arr.ptr[i]);
			transfer_info<T>::push_back(val, r);
		}
	}
};

template<class T>
struct msgpack_unwrap_impl<transfer_type_tuple, T>
{
	static void apply_tuple(const msgpack_object_array& r, size_t i, const boost::tuples::null_type& nothing) 
	{}

	template <class Head, class Tail>
	static void apply_tuple(const msgpack_object_array& r, size_t i, boost::tuples::cons<Head, Tail>& x) 
	{
		msgpack_unwrap(x.get_head(), r.ptr[i]);
		apply_tuple(r, i+1, x.get_tail());
	}

	static void unwrap(T& val, const msgpack_object& mpo)
	{
		if (mpo.type != MSGPACK_OBJECT_ARRAY)
			throw deserialization_error("msgpack object is not of type MSGPACK_OBJECT_ARRAY");
		const msgpack_object_array& arr = mpo.via.array;
		typedef typename transfer_info<T>::tuple_type tuple_type;
		if (arr.size != (size_t) boost::tuples::length<tuple_type>::value)
			throw deserialization_error("Invalid number of tuple elements");
		tuple_type out_tuple;
		apply_tuple(arr, 0, out_tuple);
		transfer_info<T>::from_boost_tuple(val, out_tuple);
	}
};

template<class T>
struct msgpack_unwrap_impl<transfer_type_map_object, T>
{
        static void unwrap(T& val, const msgpack_object& mpo)
        {
		if (mpo.type != MSGPACK_OBJECT_MAP)
			throw deserialization_error("msgpack object is not of type MSGPACK_OBJECT_MAP");

                val = T();
		const msgpack_object_map& map = mpo.via.map;
		for(size_t i = 0; i < map.size; i++)
		{
			std::string key;
			msgpack_unwrap(key, map.ptr[i].key);
			msgpack_unwrap(val[key], map.ptr[i].val);
		}
        }
};

template<class T>
void msgpack_unwrap(T& val, const msgpack_object& mpo)
{
	typedef typename transfer_info<T>::type sub_type;
	msgpack_unwrap_impl<sub_type, T>::unwrap(val, mpo);
}

template<class T, class In = std::string>
void msgpack_deserialize(T& out, const In& s)
{
	msgpack_unpacked msg;
	msgpack_unpacked_init(&msg);
	msgpack_unpack_next(&msg, s.data(), s.size(), NULL);
	try
	{
		msgpack_unwrap(out, msg.data);
	}
	catch(const io_exception& io)
	{
		msgpack_unpacked_destroy(&msg);
		throw;
	}
	msgpack_unpacked_destroy(&msg);
}

template<typename T, typename In>
T msgpack_deserialize(const In& s)
{
	T unserialized_object;
	msgpack_deserialize(unserialized_object, s);
#if __GNUC__ > 8
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wpessimizing-move"
#endif
	return std::move(unserialized_object);
#if __GNUC__ > 8
# pragma GCC diagnostic pop
#endif
}

inline void check_type(const msgpack_object& obj, int& type)
{
	if (obj.type != MSGPACK_OBJECT_ARRAY)
		throw deserialization_error("msgpack object is not of type MSGPACK_OBJECT_ARRAY");
	const msgpack_object_array& arr = obj.via.array;
	if (arr.size != 2) 
		throw deserialization_error("Type message not composed of two parts");
	if (arr.ptr[0].type != MSGPACK_OBJECT_POSITIVE_INTEGER)
		throw deserialization_error("Type message doesn't start with an interger");
	type = (int) arr.ptr[0].via.u64;
}

inline void msgpack_deserialize(msgpack_unpacked& msg, int& type, const std::string& s)
{
	msgpack_unpacked_init(&msg);
	msgpack_unpack_next(&msg, s.data(), s.size(), NULL);
	try
	{
		check_type(msg.data, type);
	}
	catch(const io_exception& io)
	{
		msgpack_unpacked_destroy(&msg);
		throw;
	}
}

template<class T>
void msgpack_decode(T& out, msgpack_unpacked& msg)
{
	msgpack_unwrap(out, msg.data.via.array.ptr[1]);
	msgpack_unpacked_destroy(&msg);
}

template<class T>
void msgpack_deserialize(T& out, int type, const std::string& s)
{
	msgpack_unpacked msg;
	int actual_type;
	msgpack_deserialize(msg, actual_type, s);
	if (actual_type != type)
	{
		msgpack_unpacked_destroy(&msg);
		throw io_exception(printstring("Type mismatch in typed deserialize (%d vs %d)", actual_type, type));
	}
	msgpack_decode(out, msg);
}

template<class SomeType>
class get_type_id
{
public:
	static const int type_id = SomeType::type_id;
};

#define SET_TYPE_ID(type, id) template<> class get_type_id<type> { public: static const int type_id = id; };

template<class T>
void msgpack_deserialize_typed(T& out, const std::string& s)
{
	msgpack_deserialize(out, get_type_id<T>::type_id, s);
}

template<class T>
void msgpack_serialize_typed(const T& in)
{
	msgpack_serialize(get_type_id<T>::type_id, in);
}

#endif

