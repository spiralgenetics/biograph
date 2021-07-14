#include "modules/web/url_query.h"
#include "modules/io/io.h"
#include "modules/web/urlencode.h" 
#include "modules/io/log.h" 

#include <utility> 
#include <vector> 
#include <boost/regex.hpp> 
#include <boost/algorithm/string.hpp>

illformed_query_string::illformed_query_string(const std::string& query)
	:io_exception("invalid query string: "+ query) {}

query_variables::query_variables() {}

query_variables::query_variables(const std::map<std::string, std::string>& that)
{
	std::map<std::string,std::string>::operator=(that);
}

// make_kv_pair("name=value") => (urldecode("name"), urldecode("value")) else throw
std::pair<std::string,std::string> make_kv_pair(const std::string& kv)
{
	std::vector<std::string> key_value_pair;
	boost::split(key_value_pair, kv, boost::is_any_of("="));
	if(key_value_pair.size() != 2)
		throw kv;
	return std::make_pair(urldecode(key_value_pair[0]) , urldecode(key_value_pair[1]));
}

query_variables::query_variables(const std::string& query)
{
	if( query.empty() )
		return;
	std::vector<std::string> components;
	boost::split(components, query, boost::is_any_of("&"), boost::token_compress_on);
	try
	{
		std::for_each( components.begin(), components.end(), [&] (const std::string& kvp) { if(!kvp.empty()) insert(make_kv_pair(kvp)); });
	}
	catch(const std::string& e)
	{
		SPLOG("unusual query string: %s <= invalid key-pair %s", query.c_str(), e.c_str());
	}
}

query_variables::operator std::string () const
{
	std::vector<std::string> kvs;
	std::for_each( cbegin(), cend(), [&kvs] ( const value_type& pair ) { kvs.push_back(pair.first+'='+pair.second); });
	return boost::join( std::move(kvs), "&");
}

bool query_variables::operator==(const query_variables& qv2) const
{
	if( size() != qv2.size() )
		return false;
	for(auto pair : *this)
	{
		query_variables::const_iterator cit = qv2.find(pair.first);
		if( cit == qv2.cend() )
			return false;
		else if( (*cit).second.compare(pair.second) == 0)
			continue;
		else
			return false;
	}
	return true;
}
