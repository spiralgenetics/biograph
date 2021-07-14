#include "modules/pipeline/dataset_path.h"
#include "modules/pipeline/ottoman.h"
#include "modules/web/urlencode.h"
#include "modules/web/couchdb.h"
#include "modules/web/httpclient.h"
#include "modules/io/config.h"

#include <boost/algorithm/string.hpp>

dataset_path::dataset_path(std::string url, bool decode)
{
	// SPLOG("dataset_path::dataset_path> %s", url.c_str());
	if (decode) {
		url = urldecode(url);
	}

	// Remove trailing slash
	if (url.size() != 0 && url[url.size()-1] == '/') {
		url = url.substr(0, url.size()-1);
	}
	// Remove leading slash (do I need to do this with token_compress_on?)
	if (url.size() != 0 && url[0] == '/') {
		url = url.substr(1);
	}

	std::vector<std::string> parts;
	boost::split(parts, url, boost::is_any_of("/"), boost::token_compress_on);

	try {
		if (parts.size() < 2) {
			throw io_exception("Invalid dataset url: too few parts");
		}
		if (parts[0] != "api") {
			throw io_exception("Invalid dataset url: no api");
		}

		size_t rest_index;
		if (parts[1] == "reference") {
			m_reference = true;
			m_user = "";
			m_base = path(CONF_S(path_reference_base));
			m_url = "/api/reference";
			rest_index = 2;
		}
		else if (parts[1] == "users") {
			if (parts.size() < 4) {
				throw io_exception("Invalid dataset url: too few parts for user URL");
			}
			if (parts[3] != "data") {
				throw io_exception("Invalid dataset url: data is not 4th part");
			}

			m_reference = false;
			m_user = parts[2];
			m_base = path(CONF_S(path_user_base)).append(m_user);
			m_url = "/api/users/" + urlencode(parts[2]) + "/data";
			rest_index = 4;
		}
		else {
			throw io_exception("Invalid dataset url: not reference or users");
		}

		m_data_dir = m_base.append("data");
		m_meta = m_base.append("meta");
		m_parent = "";
		m_name = "";

		for (size_t i = rest_index; i < parts.size(); i++) {
			auto part = parts[i];
			if (part[0] == '.') {
				throw io_exception("Invalid directory or filename");
			}
			m_parent = m_url;
			m_name = part;
			m_meta = m_meta.append(part);
			m_url += "/" + urlencode(part);
		}

		// SPLOG("dataset_path::dataset_path> parent: %s", m_parent.c_str());
		// SPLOG("dataset_path::dataset_path> name: %s", m_name.c_str());
		// SPLOG("dataset_path::dataset_path> meta: %s", m_meta.url().c_str());
		// SPLOG("dataset_path::dataset_path> url: %s", m_url.c_str());
	}
	catch (const io_exception& e) {
		throw io_exception(printstring("%s in %s", e.message().c_str(), url.c_str()));
	}
}

std::string dataset_path::friendly() const
{
	if (m_reference) {
		return CONF_S(reference_path) + "/" + m_name;
	}
	return std::string("/") + m_name;
}

dataset_path dataset_path::append(const std::string& name) const
{
	//SPLOG("dataset_path appending name=%s to for m_url=%s", name.c_str(), m_url.c_str());
	if (name == "") {
		throw io_exception("Trying to append empty name to path");
	}
	dataset_path r;
	r.m_reference = m_reference;
	r.m_url = m_url + "/" + urlencode(name);
	r.m_user = m_user;
	r.m_parent = m_url;
	//SPLOG("dataset_path url:%s appending name=%s", m_url.c_str(), m_name.c_str());
	r.m_name = name;
	r.m_base = m_base;
	r.m_meta = m_meta.append(name);
	r.m_data_dir = m_data_dir;
	return r;
}

dataset_path dataset_path::root() const
{
	dataset_path r;
	r.m_reference = m_reference;
	r.m_user = m_user;
	if (m_user != "") {
		r.m_url = "/api/users/" + urlencode(m_user) + "/data";
	}
	else {
		r.m_url = "/api/reference";
	}
	r.m_parent = "";
	r.m_name = "";
	r.m_base = m_base;
	r.m_meta = m_base.append("meta");
	r.m_data_dir = m_base.append("data");
	return r;
}

dataset_path dataset_path::root(const std::string& user)
{
	dataset_path r;
	r.m_reference = false;
	r.m_user = user;
	r.m_url = "/api/users/" + urlencode(r.m_user) + "/data";
	r.m_parent = "";
	r.m_name = "";
	r.m_base = path(CONF_S(path_user_base)).append(r.m_user);
	r.m_meta = r.m_base.append("meta");
	r.m_data_dir = r.m_base.append("data");
	return r;
}

path::exist_enum dataset_path::exists() const
{
	couch_server<direntry> db(ottoman_url());
	direntry row;

	if (db.get(row, m_url)) {
		if (row.directory) {
			return path::e_directory;
		}
		else {
			return path::e_file;
		}
	}
	return path::e_no_exist;
}

direntry dataset_path::stat() const
{
	couch_server<direntry> db(ottoman_url());
	direntry row;
	if (!db.get(row, m_url)) {
		throw io_exception(printstring("Unable to find entry for file: %s", m_url.c_str()));
	}
	return row;
}

void create_ancestors(const dataset_path& p)
{
	if (!p.parent().empty()) {
		// SPLOG("do mkdir for %s", p.parent().c_str());
		dataset_path parent_dir(p.parent().c_str());
		switch (parent_dir.exists()) {
		case path::e_no_exist:
			if (!parent_dir.parent().empty()) {
				parent_dir.mkdir();
			}
			break;
		case path::e_directory:
			break;
		case path::e_file:
			throw io_exception(
				printstring("cannot mkdir '%s' because parent path '%s' collides with existing file",
					p.friendly().c_str(), parent_dir.friendly().c_str()));
		}
	}
}

void dataset_path::create(const dataset_meta& meta) const
{
	couch_server<direntry> db(ottoman_url());
	direntry row;
	if (db.get(row, m_url)) {
		throw io_exception(printstring("Path %s already exists in dataset_path::create", m_url.c_str()));
	}
	m_meta.json_put(meta);
	set_file_direntry(*this, meta, time(0), row);
	if (!db.put(m_url, row)) {
		throw io_exception(printstring("Failed to create metadata for %s", m_url.c_str()));
	}
	create_ancestors(*this);
}

void dataset_path::update(const dataset_meta& meta) const
{
	couch_server<direntry> db(ottoman_url());
	direntry row;
	if (!db.get(row, m_url)) {
		throw io_exception(printstring("Path %s doesn't exist in dataset_path::update", m_url.c_str()));
	}
	m_meta.json_put(meta);
	set_file_direntry(*this, meta, time(0), row);
	if (!db.put(row._id, row)) {
		throw io_exception(printstring("Failed to update metadata for %s", m_url.c_str()));
	}
}

void dataset_path::remove(bool recursive) const
{
	if (recursive) {
		for (const auto& child : list_dir()) {
			dataset_path child_path(child.url);
			child_path.remove(true);
		}
	}

	couch_server<direntry> db(ottoman_url());
	direntry row;
	if (!db.get(row, m_url)) {
		SPLOG_P(LOG_DEBUG, "dataset URL not found: %s", m_url.c_str());
		return;
	}
	if (!db.erase(row)) {
		throw io_exception("Failed to do db erase");
	}

	if (row.directory) {
		m_meta.rmdir();
	}
	else {
		m_meta.remove();
	}
}

direntry dataset_path::create_remote(const dataset_meta& meta) const
{
	direntry de;
	m_meta.json_put(meta);
	set_file_direntry(*this, meta, time(0), de);
	return de;
}

void dataset_path::update_cache(const direntry& de)
{
	couch_server<direntry> db(ottoman_url());
	direntry cur;
	if (!db.get(cur, de.url)) {
		throw io_exception(printstring("Attempting to update invalid entry: %s", de.url.c_str()));
	}
	direntry next = de;
	next._id = cur._id;
	next._rev = cur._rev;
	if (!db.put(next._id, next)) {
		throw io_exception(printstring("Failed to update cache entry: %s", de.url.c_str()));
	}
}

void dataset_path::mkdir() const
{
	couch_server<direntry> db(ottoman_url());
	if (m_name.empty()) {
		throw io_exception("Empty name string in mkdir");
	}
	m_meta.mkdir();
	direntry de;
	set_directory_direntry(*this, time(0), de);
	if (!db.put(m_url, de)) {
		throw io_exception("Failed to make directory metadata");
	}
	create_ancestors(*this);
}

std::vector<direntry> dataset_path::list_dir() const
{
	couch_server<direntry> db(ottoman_url());
	return db.find_match<direntry>("by_parent", m_url);
}

void copy_dataset(const dataset_path& out, const dataset_path& in)
{
	dataset_meta dm;
	in.load(dm);
	out.create(dm);
}

void set_directory_direntry(const dataset_path& path, const time_t t, direntry& de)
{
	de.url = path.url();
	de.parent = path.parent();
	de.name = path.name();
	de.user = path.user();
	de.reference = path.is_reference();
	de.created = t;
	de.directory = true;
	de.size = 0;
	de.records = 0;
	de.in_progress = false;
}

void set_file_direntry(const dataset_path& path, const dataset_meta& meta, const time_t t, direntry& de)
{
	de.url = path.url();
	de.parent = path.parent();
	de.name = path.name();
	de.user = path.user();
	de.reference = path.is_reference();
	de.created = t;
	de.directory = false;
	de.type = meta.type;
	de.sort_keys = meta.sort_keys;
	de.size = meta.the_manifest.get_size();
	de.records = meta.the_manifest.get_num_records();
	de.ref_name = meta.ref_name;
	de.in_progress = meta.in_progress;
}

static
void create_direntry(dataset_path& path, direntry& de)
{
	time_t modtime = path.meta().modify_time();
	size_t count = 0;
	while (true) {
		try {
			dataset_meta dm;
			path.load(dm);
			set_file_direntry(path, dm, modtime, de);
			return;
		}
		catch (const io_exception& e) {
			if (count++ <= 3) {
				continue;
			}
			throw e;
		}
	}
}

static
void create_cache(http_client& couch, const dataset_path& start)
{
	std::queue<dataset_path> working;
	working.push(start);
	while (!working.empty()) {
		dataset_path path = working.front();
		working.pop();
		path::exist_enum e = path.meta().exists();
		if (e == path::e_no_exist) {
			SPLOG_P(LOG_ERR, "ERROR: %s path does not exist", path.meta().url().c_str());
			continue;
		}

		direntry de;
		if (e == path::e_directory) {
			set_directory_direntry(path, 0, de);
			auto subfiles = path.meta().list();
			for (size_t i = 0; i < subfiles.size(); i++) {
				working.push(path.append(subfiles[i]));
			}
		}
		else {
			create_direntry(path, de);
		}
		std::string key = de.url;
		std::string ignore;
		couch.do_put(urlencode(key), json_serialize(de), ignore);
	}
}

void gen_cache(const char* user)
{
	std::string couch_url = ottoman_url();
	http_client couch(couch_url);

	create_cache(couch, dataset_path("/api/reference/"));

	std::vector<std::string> users;
	if (user) {
		users.push_back(user);
	}
	else {
		path user_path(CONF_S(path_user_base));
		user_path.mkdir();
		users = user_path.list();
	}

	for (const auto& user : users) {
		create_cache(couch, dataset_path("/api/users/" + user + "/data/"));
	}
}
