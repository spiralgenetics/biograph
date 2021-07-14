#include "modules/pipeline/ottoman.h"
#include "modules/pipeline/direntry.h"
#include "modules/web/restful.h" 
#include "modules/web/couchdb.h"
#include "modules/web/httpclient.h"
#include "modules/io/log.h"

#include <mutex>

typedef couch_row<std::string, direntry> ottoman_row;
typedef couch_results<std::string, direntry> ottoman_results;
typedef std::map<std::string, direntry> ottoman_index;

const char* VIEW_URL = "/view/by_parent";

class ottoman_impl
{
public:
	ottoman_impl() {}

	bool get(const std::string& path, direntry& de);
	void put(const std::string& path, const direntry& de);
	bool del(const std::string& path);

	void by_parent(const std::string& path, ottoman_results& results);

private:
	ottoman_index m_by_path;
	std::map<std::string, ottoman_index> m_by_parent;
	std::mutex m_mutex;

	typedef std::lock_guard<std::mutex> lock_t;
};

class data_handler : public easy_rest_handler
{
public:
	data_handler(ottoman_impl& impl, http_request& request);
	std::string easy_get() override;
	bool easy_put(const std::string& input) override;
	bool easy_del() override;

private:
	std::string query();

private:
	ottoman_impl& m_impl;
};

bool ottoman_impl::get(const std::string& path, direntry& de)
{
	auto it = m_by_path.find(path);
	if (it == m_by_path.end()) {
		return false;
	}
	de = it->second;
	return true;
}

void ottoman_impl::put(const std::string& path, const direntry& de)
{
	lock_t lock(m_mutex);
	// SPLOG("ottoman_impl::put> path: %s, parent: %s", path.c_str(), de.parent.c_str());
	auto& index = m_by_parent[de.parent];
	index[path] = de;
	m_by_path[path] = de;
}

bool ottoman_impl::del(const std::string& path)
{
	lock_t lock(m_mutex);
	auto it = m_by_path.find(path);
	if (it != m_by_path.end()) {
		auto jt = m_by_parent.find(it->second.parent);
		if (jt != m_by_parent.end()) {
			auto& index = jt->second;
			index.erase(path);
			if (index.empty()) {
				m_by_parent.erase(jt);
			}
		}
		m_by_path.erase(it);
		return true;
	}
	return false;
}

void ottoman_impl::by_parent(const std::string& path, ottoman_results& results)
{
	lock_t lock(m_mutex);
	// SPLOG("ottoman_impl::by_parent> path: %s", path.c_str());
	auto it = m_by_parent.find(path);
	if (it == m_by_parent.end()) {
		// SPLOG("ottoman_impl::by_parent> index not found");
		return;
	}

	for (const auto& jt : it->second) {
		ottoman_row row;
		row.key = jt.first;
		row.value = jt.second;
		results.rows.push_back(row);
	}
	results.total_rows = results.rows.size();
}

ottoman_server::ottoman_server()
	: m_impl(new ottoman_impl)
{
	register_handler("/spiral_files(.*)", 
		[&] (http_request& request) {
			return new data_handler(*m_impl, request);
		}
	);
}

ottoman_server::~ottoman_server() {}

data_handler::data_handler(ottoman_impl& impl, http_request& request)
	: easy_rest_handler(request)
	, m_impl(impl)
{
}

std::string data_handler::easy_get()
{
	auto path = get_match_result(1);
	if (path == VIEW_URL) {
		return query();
	}

	direntry de;
	if (m_impl.get(path, de)) {
		return json_serialize(de);
	}

	return "";
}

bool data_handler::easy_put(const std::string& input)
{
	auto path = get_match_result(1);
	if (path == VIEW_URL) {
		return false;
	}

	direntry de;
	json_deserialize(de, input);
	de._id = path;
	m_impl.put(path, de);

	return true;
}

bool data_handler::easy_del()
{
	auto path = get_match_result(1);
	if (path == VIEW_URL) {
		return false;
	}

	m_impl.del(path);
	return true;
}

std::string data_handler::query()
{
	auto key = m_request.get_variable("key");
	// SPLOG("data_handler::query> key: %s", key.c_str());
	if (key.empty()) {
		return "";
	}
	std::string parent;
	json_deserialize(parent, key);
	ottoman_results results;
	m_impl.by_parent(parent, results);
	auto json = json_serialize(results);
	// SPLOG("data_handler::query> json: %s", json.c_str());
	return json;
}

std::string ottoman_url()
{
	return make_client_url(
		"ottoman_bind_list",
		"MASTER_PORT_5984_TCP_ADDR", 
		"MASTER_PORT_5984_TCP_PORT", 
		"/spiral_files/"
	);
}
