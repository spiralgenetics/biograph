#include "modules/mapred/metadata.h"

const char* meta::ns::internal = "internal";
const char* meta::ns::readonly = "spiral_readonly";
const char* meta::ns::user     = "spiral";

// All metadata so far:
//  namespace   key                        type
//  ---------   ---                        ----
//  internal    encoding                   string
//  internal    created                    time_t
//  internal    kmer_db                    string (path)
//  internal    entries                    size_t
//  readonly    read_size                  size_t
//  readonly    created                    string (RFC 3339)
//  readonly    kmer_size                  size_t
//  readonly    sample_id                  string (UUID or MD5 hex)
//  readonly    sample_bases               size_t
//  readonly    corrected_read_count       size_t
//  readonly    corrected_base_dist        vector<size_t>
//  readonly    failed_correction_count    size_t
//  readonly    filtered_kmers             size_t
//  readonly    processed_read_count       size_t

static
std::map<std::string, meta::merge::function>& registry()
{
 	static std::map<std::string, meta::merge::function> registry;
 	return registry;
}

meta::merge::init meta::merge::register_fn(const std::string& key, meta::merge::function fn)
{
	auto result = registry().insert(std::make_pair(key, fn));
	if (!result.second) {
		throw io_exception(printstring("meta::merge::register_fn> already registered: %s ", 
			key.c_str()
		));
	}
	return init{};
}

void meta::data::set_creation_time_now()
{
	auto now = std::time(nullptr);
	set(meta::ns::readonly, "created", time_to_RFC3339(now));
	set(meta::ns::internal, "created", now);
}

void meta::data::set_options(const std::string& step_name, const std::string& the_options_json)
{
	set(meta::ns::readonly, step_name + "_parameters", the_options_json);
}

void meta::data::set_runtime(time_t start_time)
{
	set(meta::ns::readonly, "wall_clock_runtime_sec", std::time(nullptr) - start_time);
}

void meta::data::merge(const meta::data& in)
{
	for (auto& that_ns : in.m_data) {
		auto ns_it = m_data.find(that_ns.first);
		if (ns_it == m_data.cend()) {
			// this metadata does not have this incoming namespace
			// let's copy the entire map
			m_data[that_ns.first] = that_ns.second;
		}
		else {
			for (auto that_kv : that_ns.second) {
				auto kv_it = ns_it->second.find(that_kv.first);
				if (kv_it == ns_it->second.cend()) {
					ns_it->second[that_kv.first] = that_kv.second;
				}
				else {
					// perform merge: same (ns, key)
					meta::merge::function merger{meta::merge::collide};
					auto fn_it = registry().find(kv_it->first);
					if (fn_it != registry().end()) {
						merger = fn_it->second;
					}
					meta::merge::params params{
						.ns = ns_it->first,
						.key = kv_it->first,
						.value1 = kv_it->second,
						.value2 = that_kv.second,
					};
					ns_it->second[kv_it->first] = merger(params);
				}
			}
		}
	}
}

void meta::data::unset(const std::string& ns, const std::string& key)
{
	auto ns_it = m_data.find(ns);
	if (ns_it != m_data.cend()) {
		auto key_it = ns_it->second.find(key);
		if (key_it != ns_it->second.cend()) {
			ns_it->second.erase(key_it);
		}
	}
}

js::mValue meta::merge::first(const params& params)
{
	return params.value1;
}

js::mValue meta::merge::second(const params& params)
{
	return params.value2;
}

js::mValue meta::merge::sum(const params& params)
{
	if (params.value1.type() != params.value2.type()) {
		throw io_exception("meta::merge::sum> mismatched types");
	}

	switch (params.value1.type()) {
	case js::obj_type:
		throw io_exception("meta::merge::sum> unsupported obj_type");
	case js::array_type:
		throw io_exception("meta::merge::sum> unsupported array_type");
	case js::str_type:
		throw io_exception("meta::merge::sum> unsupported str_type");
	case js::bool_type:
		throw io_exception("meta::merge::sum> unsupported bool_type");
	case js::int_type:
		return params.value1.get_uint64() + params.value2.get_uint64();
	case js::real_type:
		return params.value1.get_real() + params.value2.get_real();
	case js::null_type:
		throw io_exception("meta::merge::sum> unsupported null_type");
	}

	throw io_exception("meta::merge::sum> unknown type");
}

js::mValue meta::merge::collide(const params& params)
{
	if (!(params.value1 == params.value2)) {
		std::string str1 = json_serialize(params.value1);
		std::string str2 = json_serialize(params.value2);
		throw io_exception(printstring("Metadata collision detected for %s/%s: { %s != %s }",
			params.ns.c_str(),
			params.key.c_str(),
			str1.c_str(),
			str2.c_str()
		));
	}
	return params.value1;
}

meta::merge::init merge_created(meta::merge::register_fn("created", [] (
		const meta::merge::params& params)
	{
		auto now = std::time(nullptr);
		if (params.ns == meta::ns::internal) {
			return js::mValue(now);
		}
		return js::mValue(time_to_RFC3339(now));
	}
));

meta::merge::init merge_encoding(meta::merge::register_fn("encoding", meta::merge::first));
meta::merge::init merge_wall_clock_runtime_sec(meta::merge::register_fn("wall_clock_runtime_sec", meta::merge::second));
meta::merge::init merge_kmer_db(meta::merge::register_fn("kmer_db", meta::merge::second));
meta::merge::init merge_kmer_filter_params(meta::merge::register_fn("kmer_filter_parameters", meta::merge::second));
meta::merge::init merge_filtered_kmers(meta::merge::register_fn("filtered_kmers", meta::merge::second));
meta::merge::init merge_filtered_read_dist(meta::merge::register_fn("filtered_read_dist", meta::merge::second));
meta::merge::init merge_filtered_assembly_dist(meta::merge::register_fn("filtered_assembly_dist", meta::merge::second));
meta::merge::init merge_tagged_read_count(meta::merge::register_fn("tagged_reads_count", meta::merge::second));
meta::merge::init merge_tagged_assembly_count(meta::merge::register_fn("tagged_assembly_count", meta::merge::second));
