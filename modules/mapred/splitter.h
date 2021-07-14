#ifndef __MAPRED_SPLITTER_H__
#define __MAPRED_SPLITTER_H__

#include <string>
#include "modules/io/registry.h"

// Used by the chunker to determine whether a new file should be split off.
class splitter
{
public:
	virtual ~splitter();
	
	// Returns true when a new file should be split.
	virtual bool operator() (const std::string& key) const = 0;
	virtual void set_initial_key(const std::string& key) {}
};


// A placeholder splitter that never splits.
class null_splitter : public splitter
{
public:
	null_splitter(const std::string& /*key*/) {}
	virtual ~null_splitter() {}
	bool operator() (const std::string& key) const override { return false; }
	void set_initial_key(const std::string& key) override {}
};

DECLARE_REGISTRY_1(splitter, std::string const&);

#endif // __MAPRED_SPLITTER_H__
