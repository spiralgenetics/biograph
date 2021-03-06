#pragma once

#include "modules/io/io.h"
#include "modules/io/utils.h" 

#include <stdint.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/container/small_vector.hpp>

class deserialization_error : public io_exception
{
public:
	deserialization_error(const std::string& err) 
		: io_exception(printstring("deserialization: %s",err.c_str()))
	{}
};

class transfer_flags
{
public:
	explicit transfer_flags(int flag_val) 
	{ 
		m_val = flag_val; 
	}

	transfer_flags operator|(const transfer_flags& rhs) const
	{ 
		return transfer_flags(m_val | rhs.m_val); 
	}

	bool operator&(const transfer_flags& rhs) const 
	{ 
		return (m_val & rhs.m_val); 
	}

private:
	int m_val;
};

#define TF_STRICT	transfer_flags(1)   /* Complain if field is missing */   
#define TF_ALLOW_NULL	transfer_flags(2)   /* Don't complain if null, just apply default */
#define TF_NO_DEFAULT	transfer_flags(4)   /* Don't default construct missing fields */

typedef int64_t transfer_type_signed;
typedef uint64_t transfer_type_unsigned;
typedef double transfer_type_real;
typedef bool transfer_type_boolean;
typedef std::string transfer_type_string;
struct transfer_type_object {};
struct transfer_type_map_object {};
struct transfer_type_array {};
struct transfer_type_tuple {};
struct transfer_type_null {};

template <class Obj>
struct transfer_info
{
	typedef transfer_type_object type;
	template <class Context>
	static void object_transfer(Context& ctx, Obj& obj)
	{
		obj.transfer(ctx);
	}
};

// Automatically 'deconst' objects
template <class Obj>
struct transfer_info<const Obj>
{
	typedef transfer_type_object type;
	template <class Context>
	static void object_transfer(Context& ctx, const Obj& obj)
	{
		const_cast<Obj&>(obj).transfer(ctx);
	}
};

#define BASE_TYPE(native_type, transfer_type) \
template <>  \
struct transfer_info<native_type> \
{ \
	typedef transfer_type type; \
	static transfer_type get(const native_type& x) { return (transfer_type) x; } \
	static void put(native_type& x, const transfer_type& tt) { x = (native_type) tt; } \
};

BASE_TYPE(char, transfer_type_signed);
BASE_TYPE(unsigned char, transfer_type_unsigned);
BASE_TYPE(short, transfer_type_signed);
BASE_TYPE(unsigned short, transfer_type_unsigned);
BASE_TYPE(int, transfer_type_signed);
BASE_TYPE(unsigned int, transfer_type_unsigned);
BASE_TYPE(int64_t, transfer_type_signed);
BASE_TYPE(uint64_t, transfer_type_unsigned);
//BASE_TYPE(size_t, transfer_type_unsigned);
BASE_TYPE(bool, transfer_type_boolean);
BASE_TYPE(float, transfer_type_real);
BASE_TYPE(double, transfer_type_real);
BASE_TYPE(std::string, transfer_type_string);
BASE_TYPE(transfer_type_null, transfer_type_null);

template <class Inner>
struct transfer_info<std::vector<Inner>>
{
	typedef transfer_type_array type;
	typedef Inner value_type;
	typedef typename std::vector<Inner>::iterator iterator;
	static size_t size(std::vector<Inner>& arr) { return arr.size(); }
	static iterator begin(std::vector<Inner>& arr) { return arr.begin(); }
	static iterator end(std::vector<Inner>& arr) { return arr.end(); }
	static void push_back(std::vector<Inner>& arr, const Inner& inner) { arr.push_back(inner); }
};

template <class Inner, size_t N>
struct transfer_info<boost::container::small_vector<Inner, N>>
{
	typedef transfer_type_array type;
	typedef Inner value_type;
	typedef typename boost::container::small_vector<Inner, N>::iterator iterator;
	static size_t size(boost::container::small_vector<Inner, N>& arr) { return arr.size(); }
	static iterator begin(boost::container::small_vector<Inner, N>& arr) { return arr.begin(); }
	static iterator end(boost::container::small_vector<Inner, N>& arr) { return arr.end(); }
	static void push_back(boost::container::small_vector<Inner, N>& arr, const Inner& inner) { arr.push_back(inner); }
};

template <class Inner>
struct transfer_info<std::set<Inner>>
{
	typedef transfer_type_array type;
	typedef Inner value_type;
	typedef typename std::set<Inner>::iterator iterator;
	static size_t size(std::set<Inner>& arr) { return arr.size(); }
	static iterator begin(std::set<Inner>& arr) { return arr.begin(); }
	static iterator end(std::set<Inner>& arr) { return arr.end(); }
	static void push_back(std::set<Inner>& arr, const Inner& inner) { arr.insert(inner); }
};

template <class T1, class T2>
struct transfer_info<std::pair<T1, T2>>
{
	typedef transfer_type_tuple type;
	typedef boost::tuple<T1, T2> tuple_type;
	static tuple_type as_boost_tuple(const std::pair<T1, T2>& pair) 
	{ return boost::make_tuple(pair.first, pair.second); }
	static void from_boost_tuple(std::pair<T1, T2>& pair, const tuple_type& t) 
	{ pair = std::make_pair(boost::get<0>(t), boost::get<1>(t)); }
};

// Special version of pair to 'de-const' first element 
// Used in map
template <class T1, class T2>
struct transfer_info<std::pair<const T1, T2>>
{
	typedef transfer_type_tuple type;
	typedef boost::tuple<T1, T2> tuple_type;
	static tuple_type as_boost_tuple(const std::pair<const T1, T2>& pair) 
	{ return boost::make_tuple(pair.first, pair.second); }
	static void from_boost_tuple(std::pair<const T1, T2>& pair, const tuple_type& t) 
	{ pair = std::make_pair(boost::get<0>(t), boost::get<1>(t)); }
};

template <class K, class V>
struct transfer_info<std::map<K, V>>
{
	typedef transfer_type_array type;
	typedef std::pair<K, V> value_type;
	typedef typename std::map<K,V>::iterator iterator;
	static size_t size(std::map<K, V>& m) { return m.size(); }
	static iterator begin(std::map<K, V>& m) { return m.begin(); }
	static iterator end(std::map<K, V>& m) { return m.end(); }
	static void push_back(std::map<K, V>& m, const value_type& inner) { m.insert(inner); }
};

template <class K, class V>
struct transfer_info<std::unordered_map<K, V>>
{
	typedef transfer_type_array type;
	typedef std::pair<K, V> value_type;
	typedef typename std::unordered_map<K, V>::iterator iterator;
	static size_t size(std::unordered_map<K, V>& m) { return m.size(); }
	static iterator begin(std::unordered_map<K, V>& m) { return m.begin(); }
	static iterator end(std::unordered_map<K, V>& m) { return m.end(); }
	static void push_back(std::unordered_map<K, V>& m, const value_type& inner) { m.insert(inner); }
};

template <class K, class V, class H>
struct transfer_info<std::unordered_map<K, V, H>>
{
	typedef transfer_type_array type;
	typedef std::pair<K, V> value_type;
	typedef typename std::unordered_map<K, V, H>::iterator iterator;
	static size_t size(std::unordered_map<K, V, H>& m) { return m.size(); }
	static iterator begin(std::unordered_map<K, V, H>& m) { return m.begin(); }
	static iterator end(std::unordered_map<K, V, H>& m) { return m.end(); }
	static void push_back(std::unordered_map<K, V, H>& m, const value_type& inner) { m.insert(inner); }
};

template <class V>
struct transfer_info<std::map<std::string, V>>
{
	typedef transfer_type_map_object type;
	typedef V value_type;
};

#define TF_STRICT	transfer_flags(1)   /* Complain if field is missing */   
#define TF_ALLOW_NULL	transfer_flags(2)   /* Don't complain if null, just apply default */
#define TF_NO_DEFAULT	transfer_flags(4)   /* Don't default construct missing fields */

template <class Context, class Object>
void transfer_field(
	Context& ctx, 
	const std::string& name, 
	int tag, 
	Object& obj, 
	const Object& def, 
	const transfer_flags& flags)
{
	if (ctx.is_serialize()) {
		ctx.transfer_field(name, tag, obj);
	}
	else {
		if (!ctx.has_field(name, tag)) {
			if (flags & TF_STRICT) {
				throw deserialization_error(printstring("Field '%s' is missing and strict is set", name.c_str()));
			}
			if (flags & TF_NO_DEFAULT) {
				return;
			}
			obj = def;
		}
		else {
			if ((flags & TF_ALLOW_NULL) && ctx.is_null(name, tag)) {
				if (!(flags & TF_NO_DEFAULT)) {
					obj = def;
				}
				return;
			}
			ctx.transfer_field(name, tag, obj);
		}
	}
}

template <class Context, class Object>
void transfer_field(Context& ctx, const std::string& name, int tag, Object& obj, const Object& def)
{
	transfer_field(ctx, name, tag, obj, def, transfer_flags(0));
}

template <class Context, class Object>
void transfer_field(Context& ctx, const std::string& name, int tag, Object& obj, const transfer_flags& flags)
{
	transfer_field(ctx, name, tag, obj, Object(), flags);
}

template <class Context, class Object>
void transfer_field(Context& ctx, const std::string& name, int tag, Object& obj)
{
	transfer_field(ctx, name, tag, obj, Object(), transfer_flags(0));
}


#define TRANSFER_OBJECT \
	template <class Context> \
	void transfer(Context& _ctx)

#define VERSION(v) int _next_tag __attribute__((unused)) = 1; _ctx.set_version(v)
#define FIELD(name, ...) transfer_field(_ctx, #name, _next_tag++, name, ##__VA_ARGS__)
#define FIELD_SPECIAL(special_type, name, ...) special_type(name, ##__VA_ARGS__) 
#define OBSOLETE_FIELD(name, ...) _next_tag++
#define IS_SERIALIZE (_ctx.is_serialize())
#define IS_DESERIALIZE (!_ctx.is_serialize())
#define IS_HUMAN_READABLE (_ctx.is_human_readable())
#define GET_VERSION(v) (_ctx.get_version())

// Special types
#define COUCHDB_RESERVED(name) do { \
	if (!_ctx.is_serialize() || name != "") \
		transfer_field(_ctx, #name, _next_tag, name); \
	_next_tag++; \
	} while(0)
