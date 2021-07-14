#pragma once

#include "modules/bio_format/exporter.h"

#include <math.h>

class stats
{
public:
	stats() 
		: m_min(0)
		, m_max(0)
		, m_count(0)
		, m_total(0)
		, m_total_sq(0) 
	{}

	void add(uint64_t val, uint64_t count)
	{
		m_data[val] += count;
		if (m_total == 0 || val < m_min)
			m_min = val;
		if (m_total == 0 || val > m_max)
			m_max = val;
		m_count += count;
		m_total += count * val;
		m_total_sq += count * val * val;
	}

	uint64_t count() const { return m_count; }
	double average() const { if (m_count == 0) return 0.0; return double(m_total) / double(m_count); }
	double std_dev() const
	{
		if (m_count == 0) return 0.0;
		return sqrt(double(m_count) * double(m_total_sq) - pow(double(m_total), 2)) / (double)m_count;
	}
	uint64_t min() const { return m_min; }
	uint64_t max() const { return m_max; }
	double xtile(double perc) const {
		uint64_t tot_count = 0;
		for(auto kvp : m_data) {
			if (double(tot_count + kvp.second) / double(m_count) >= perc) {
				return kvp.first;
			}
			tot_count += kvp.second;
		}
		return -1.0;
	}
private:
	uint64_t m_min;
	uint64_t m_max;
	uint64_t m_count;
	uint64_t m_total;
	uint64_t m_total_sq;
	std::map<uint64_t, uint64_t> m_data;
};

inline void print_stats(const stats& s, writable& w)
{
	w.print("{ \"cnt\": %ld, \"avg\":%f, \"std\":%f, \"p05\": %f, \"p25\": %f, \"p50\": %f, \"p75\": %f, \"p95\": %f }",
		s.count(), s.average(), s.std_dev(), s.xtile(.05), s.xtile(.25), s.xtile(.5), s.xtile(.75), s.xtile(.95));
	/*
	w.print("%s: count=%lld, avg=%f, std-dev=%f, min=%lld, max=%lld\n", 
		prefix.c_str(), s.count(), s.average(), s.std_dev(), s.min(), s.max());
	*/
}

class writable;
class read_qual_exporter : public exporter
{
public:
	read_qual_exporter(writable & byte_sink)
		: exporter(byte_sink)
		, m_by_base(6)
		, m_by_pos(256)
		, m_by_qual(256) 
	{}
	read_qual_exporter(writable & byte_sink, bool /*unused*/, const std::string& /*serialized_data*/)
		: exporter(byte_sink)
		, m_by_base(6)
		, m_by_pos(256)
		, m_by_qual(256) 
	{}

	void write(const std::string& key, const std::string& value) override;
	void write_header() override;
	void write_footer() override;


private:
	std::vector<stats> m_by_base;
	std::vector<stats> m_by_pos;
	std::vector<stats> m_by_qual;
	stats m_all;
	stats m_ends;
};
