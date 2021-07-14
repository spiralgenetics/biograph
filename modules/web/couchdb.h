#pragma once

#include "modules/web/httpclient.h"
#include "modules/web/urlencode.h"
#include "modules/io/json_transfer.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"

/* JSON row returned by a couchdb view, we ignore the 'include_docs' option */
template<class Key, class Value>
struct couch_row
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(key, TF_STRICT);
		FIELD(value, TF_STRICT);
	}

	Key key;  // The 'key' for this row, basically the result for the 'map' function
	Value value;  // The value of the map, in this case, just a revision number
};

template<class Key, class Value>
struct couch_results
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(total_rows, TF_STRICT);
		FIELD(offset, TF_STRICT);
		FIELD(rows, TF_STRICT);
	}

	size_t total_rows = 0;  // Metaata returned by couch
	size_t offset = 0;  // Metaata returned by couch
	std::vector<couch_row<Key, Value>> rows;  // Actual row data from query
};

template<class Key>
class couch_query
{
public:
	couch_query(const std::string& index)
		: m_index(index)
	{}

	void set_begin_obj(const std::string& obj)
	{
		add_param("startkey_docid", obj);
	}

	void set_limit(size_t limit)
	{
		add_param("limit", printstring("%d", (int) limit));
	}

	void set_descending()
	{
		add_param("descending", "true");
	}

	void set_group()
	{
		add_param("group", "true");
	}

	void set_group_level(int level)
	{
		add_param("group_level", printstring("%d", level));
	}

	const std::string& get_index() const { return m_index; }
	const std::string& get_query_string() const { return m_query_string; }

	void set_key(const Key& k)
	{
		this->add_param("key", json_serialize(k));
	}

	void set_begin_key(const Key& k)
	{
		this->add_param("startkey", json_serialize(k));
	}

	void set_end_key(const Key& k)
	{
		this->add_param("endkey", json_serialize(k));
		this->add_param("inclusive_end", "false");
	}

private:
	void add_param(const std::string& key, const std::string& value)
	{
		m_query_string += (m_first ? '?' : '&');
		m_query_string += key + "=" + urlencode(value);
		m_first = false;
	}

private:
	std::string m_index;
	std::string m_query_string;
	bool m_first = true;
};

template<class Doc>
class couch_server
{
public:
	// db_url is the full URL to the database (ie, http://localhost:5984/some_db)
	couch_server(const std::string& db_url)
		: m_server(db_url)
	{}

	template<class Value, class Key>
	std::vector<Value> run_query(const couch_query<Key>& query) const
	{
		std::string url = printstring("view/%s%s",
			query.get_index().c_str(), query.get_query_string().c_str()
		);
		std::string string_result;
		couch_results<Key, Value> results;
		int status = m_server.do_get(url, string_result);
		if (status != 200) {
			throw io_exception("Invalid response from taskdb: " + string_result);
		}
		json_deserialize(results, string_result);

		std::vector<Value> out;
		out.reserve(results.rows.size());
		for (const auto& row : results.rows) {
			out.push_back(row.value);
		}

		return out;
	}

	// Helpers
	// Finds records in a view that are an exact match for a key
	template<class Value, class Key>
	std::vector<Value> find_match(const std::string& index, const Key& key, size_t limit = 0) const
	{
		couch_query<Key> query(index);
		query.set_key(key);
		if (limit) {
			query.set_limit(limit);
		}
		return run_query<Value>(query);
	}

	// Finds all values that haves keys within a specified range
	// If group = true, limit reduction to elements with the same key
	template<class Value, class Key>
	std::vector<Value> find_range(
		const std::string& index,
		const Key& start,
		const Key& end,
		size_t limit = 0,
		int group_level = 0) const
	{
		couch_query<Key> query(index);
		query.set_begin_key(start);
		query.set_end_key(end);
		if (limit) query.set_limit(limit);
		if (group_level != 0) {
			if (group_level < 0) {
				query.set_group();
			}
			else {
				query.set_group_level(group_level);
			}
		}
		return run_query<Value>(query);
	}

	// Get a document by master index, false if no such doc
	bool get(Doc& out, const std::string& key) const
	{
		// SPLOG_P(LOG_DEBUG, "couch_server::get> about to get %s", key.c_str());
		std::string tmp;
		int result = m_server.do_get(urlencode(key), tmp);
		// SPLOG_P(LOG_DEBUG, "couch_server::get> got status: %s", std::to_string(result).c_str());
		// SPLOG_P(LOG_DEBUG, "couch_server::get> got body: %s", tmp.c_str());
		if (result == 204 || result == 404) {
			return false;
		}
		if (result != 200) {
			throw io_exception(printstring("Unexpected http status code %d during taskdb get: %s", result, tmp.c_str()));
		}
		json_deserialize(out, tmp);
		return true;
	}

	// Get a document by master index, throw if no such doc
	Doc get(const std::string& key) const
	{
		Doc doc;
		if (!get(doc, key)) {
			throw io_exception(printstring("Unknown record in taskdb get for key: %s", key.c_str()));
		}
		return doc;
	}

	// Inserts a new document, return false if conflict
	bool put(const std::string& key, const Doc& val) const
	{
		//SPLOG("Inserting doc: %s at %s", json_serialize(val).c_str(), key.c_str());
		std::string ignore;
		int result = m_server.do_put(urlencode(key), json_serialize(val), ignore);
		return check_status(result, 201);
	}

	// Trys to update a document, returns false if conflict, throws on error
	bool put(const Doc& doc) const
	{
		//SPLOG("Updating doc: %s", json_serialize(doc).c_str());
		std::string ignore;
		int result = m_server.do_put(urlencode(doc._id), json_serialize(doc), ignore);
		return check_status(result, 201);
	}

	// Trys to erase a document, returns false if conflict, throws on error
	bool erase(const Doc& doc) const
	{
		std::string del_url = urlencode(doc._id) + "?rev=" + doc._rev;
		//SPLOG("Deleteing %s", del_url.c_str());
		int result = m_server.do_delete(urlencode(doc._id) + "?rev=" + doc._rev);
		return check_status(result, 200);
	}

private:
	bool check_status(int result, int goal_status) const
	{
		if (result == 409 || result == 404) {
			return false;
		}
		if (result == goal_status) {
			return true;
		}
		throw io_exception(printstring("Unexpected http status code %d during taskdb operation", result ));
	}

	mutable json_client m_server;
};
