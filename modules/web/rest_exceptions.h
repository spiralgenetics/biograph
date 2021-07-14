#pragma once

class rest_exception : public io_exception // TODO: Rename io_exception to spiral_exception
{
public:
	rest_exception(const std::string& message, int errcode)
		: io_exception(message)
		, m_errcode(errcode)
	{}

	int get_errcode() const { return m_errcode; }
private:
	int m_errcode;
};

class bad_request : public rest_exception
{
public:
	bad_request(const std::string& uri, const std::string& data="") 
		: rest_exception("Bad request " + uri + " with data: "+data, 400)
	{}
};

class method_not_allowed : public rest_exception
{
public:
	method_not_allowed(const std::string method) 
		: rest_exception("method " + method + " not allowed on this URI", 405)
	{}
};

// TODO: I should not be returning 'forbidden' for unathorized
class forbidden : public rest_exception
{
public:
	forbidden(const std::string& message)
		: rest_exception(message, 403) 
	{}
};

class uri_not_found : public rest_exception
{
public:
	uri_not_found(const std::string& uri) 
		: rest_exception("URI " + uri + " not found", 404)
	{}
};

class conflict : public rest_exception
{
public:
	conflict() : rest_exception("Conflict", 409) {}
};

class unprocessable_entity : public rest_exception
{
public:
	unprocessable_entity(const std::string& message)
		: rest_exception(message, 422) 
	{}
};

class locked : public rest_exception
{
public:
	locked(const std::string& message)
		: rest_exception(message, 423) 
	{}
};
