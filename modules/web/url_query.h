#ifndef __url_query_h__
#define __url_query_h__

#include "modules/io/io.h" 
#include <map>
#include <string>

class illformed_query_string : public io_exception
{
	public:
		illformed_query_string(const std::string& query);
};

class query_variables : public std::map<std::string,std::string>
{
	public:
		query_variables();
		query_variables(const std::map<std::string, std::string>&);
		query_variables(const std::string&);	// converts "name=value&foo=bar" into { {"name", "value"}, {"foo", "bar"} }

		operator std::string () const;		// converts { {"name", "value"}, {"foo", "bar"} } into "name=value&foo=bar"

		bool operator==(const query_variables& qv2) const;
};

#endif
