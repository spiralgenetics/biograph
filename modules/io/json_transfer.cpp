#include "modules/io/json_transfer.h"

static const std::map<js::Value_type, std::string> g_type_to_str {
	{js::obj_type,	"object"},
	{js::array_type,"array"},
	{js::str_type,	"string"},
	{js::bool_type,	"bool"},
	{js::int_type,	"int"},
	{js::real_type,	"real"},
	{js::null_type,	"null"}
};

std::string exception_msg(const js::Value_type& t)
{
	return printstring("unknown json_spirit type with enum %d", t);
}

void throw_bad_type(const js::Value_type& found_type, const js::Value_type& expected_type)
{
	auto found_it = g_type_to_str.find(found_type);
	auto expected_it = g_type_to_str.find(expected_type);
	if (found_it == g_type_to_str.end()) {
		throw deserialization_error(exception_msg(found_type));
	}
	if (expected_it == g_type_to_str.end()) {
		throw deserialization_error(exception_msg(expected_type));
	}

	std::string found_type_as_str = found_it->second;
	std::string expected_type_as_str = expected_it->second;
	throw deserialization_error("Json is of type " + found_type_as_str + ", not an " + expected_type_as_str);
}
