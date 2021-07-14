#pragma once

#include "modules/web/restful.h"
#include "modules/io/json_transfer.h"
#include "modules/io/utils.h"

#include <map>

// This class represents a 'registry' for a staticly defined set of elements
// of type D.  It also automatically provides a restful data handler to let the outside
// world see the registered elements.  It also provides a mechanism for finding elements
// by id, and a 'reference' type that holds references and writes them via URLs
// The reference type also manages static initalization order....
// D must implement an id member, a URL member, and be JSON serializable

template<class D>
class restful_registry
{
public:
	typedef std::shared_ptr<D> ref_type;
	typedef std::map<std::string, ref_type> map_t;
	
	class base_rest_handler : public easy_rest_handler
	{
	public:
		base_rest_handler(http_request& request) 
			: easy_rest_handler(request) 
		{}

		std::string easy_get() override
		{
			std::vector<js::mValue> data;
			for (const auto& item : s_map) {
				data.push_back(json_wrap(*item.second));
			}
			return json_serialize(data);
		}
	};

	static 
	rest_handler* make_base(http_request& request) { return new base_rest_handler(request); }

	class single_rest_handler : public easy_rest_handler
	{
	public:
		single_rest_handler(http_request& request) 
			: easy_rest_handler(request) 
		{}

		std::string easy_get() override
		{
			std::string key = get_match_result(1);
			auto it = restful_registry::s_map.find(key);
			if (it == restful_registry::s_map.end()) {
				throw uri_not_found(key);
			}
			return json_serialize(*it->second);
		}
	};	

	static 
	rest_handler* make_single(http_request& request) { return new single_rest_handler(request); }

	static 
	void add(std::shared_ptr<D> x)
	{
		x->url = s_base_url + "/" + x->id;
		s_map.insert(std::make_pair(x->id, x));
	}

	static 
	ref_type find(const std::string& id)
	{
		auto it = s_map.find(id);
		if (it != s_map.end()) {
			return ref_type(it->second);
		}
		throw io_exception(std::string("Cannot find type ")+id+std::string(" in the registry"));
	}

	static 
	void rest_register(const std::string& url)
	{
		s_base_url = url;
		register_handler(url, make_base);
		register_handler(url + "/(.*)", make_single);
	}

	static const std::string& base_url() { return s_base_url; }

private:
	friend class base_rest_handler;
	friend class single_rest_handler;

	static std::string s_base_url;
	static map_t s_map;
};

template<class D>
typename restful_registry<D>::map_t restful_registry<D>::s_map;

template<class D>
std::string restful_registry<D>::s_base_url;

template<class D>
struct transfer_info<std::shared_ptr<D>>
{
	typedef std::string type;
	static std::string get(const std::shared_ptr<D>& value)
	{ 
		if (value) {
			return value->url; 
		}
		return "";
	}

	static 
	void put(std::shared_ptr<D>& value, const std::string& url)
	{
		if (url == "") {
			value = std::shared_ptr<D>();
			return;
		}

		std::string base_url = restful_registry<D>::base_url();
		if (url.substr(0, base_url.size() + 1) != base_url + "/") {
			throw io_exception(printstring("Unknown URL reference: %s", url.c_str()));
		}
		std::string id = url.substr(base_url.size() + 1, url.size() - (base_url.size() + 1));
		value = restful_registry<D>::find(id);
		if (!value) {
			throw io_exception(printstring("Unknown id component: %s, in URL: %s ", id.c_str(), url.c_str()));
		}
	}
};
