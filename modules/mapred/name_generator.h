#pragma once

#include "modules/mapred/path.h"

class name_generator
{
public:
	virtual ~name_generator() {}
	virtual path generate_name(const std::string& unique) const = 0;
};

class simple_name_generator : public name_generator
{
public:
	inline simple_name_generator(const path& p) : m_p(p) {}
	inline path generate_name(const std::string& unique) const override { return m_p.append(unique); }
private:
	path m_p;
};

class prefix_name_generator : public name_generator
{
public:
	inline prefix_name_generator(const name_generator& fng, const std::string& prefix)
		: m_fng(fng)
		, m_prefix(prefix)
	{}
	inline path generate_name(const std::string& unique) const override
	{
		return m_fng.generate_name(m_prefix + unique);
	}
private:
	const name_generator& m_fng;
	std::string m_prefix;
};
