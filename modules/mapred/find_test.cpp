#include "gtest/gtest.h"
#include "modules/test/test_utils.h"
#include "modules/mapred/path.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/query.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/map_task.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"

typedef std::map<std::string, std::string> StringMap;

void genSortedData( manifest& sortedData, StringMap& verify, const size_t totalCount)
{
	//create a kv file full of random keys and values
	manifest randomData;
	path test_path(make_path("find"));
	gen_random_kv(test_path, totalCount, 1024*1024, 10, verify, randomData, codec::null);
	size_t chunk_size = 500; // Room for about 10 values

	//sort the file
	task_mgr_local tm;
	std::unique_ptr<map_task> mt = make_unique<map_task>();
	mt->input = randomData;
	mt->map = "identity";
	mt->output_goal_size = chunk_size;
	mt->sort = "lexical";//this tells the mapper to perform a sort after the map operation 'identity'

	manifest map_manifest;
	tm.run_task(map_manifest, test_path.append("do_map"), std::move(mt));
	// map_manifest now contains data locally sorted in its chunks.

	std::unique_ptr<sort_task> st = make_unique<sort_task>();
	st->input = map_manifest;
	st->goal_size = chunk_size;
	st->max_files = 8;

	// sort_task will produce a globally sorted manifest
	tm.run_task(sortedData, test_path.append("do_sort"), std::move(st));
}

TEST(find, manifest)
{
	StringMap verify;
	manifest sortedData;
	size_t totalCount = 100;
	genSortedData(sortedData, verify, totalCount);
	ASSERT_EQ( totalCount, sortedData.get_num_records() ); 

	query qr;
	std::string key,value;

	// check basic stuff
	qr.find(sortedData, "", "");
	ASSERT_FALSE(qr.read(key,value));

	std::string evenSmallerKey, smallerKey, minKey, maxKey, biggerKey, evenBiggerKey;
	maxKey = std::max_element(verify.begin(), verify.end(), verify.value_comp())->first;
	minKey = std::min_element(verify.begin(), verify.end(), verify.value_comp())->first;
	smallerKey = minKey;
	smallerKey[smallerKey.size() -1] = --smallerKey[smallerKey.size() - 1];
	evenSmallerKey = smallerKey;
	evenSmallerKey[smallerKey.size() -1] = --evenSmallerKey[smallerKey.size() - 1];
	biggerKey = maxKey+"A"; //guaranteed to be outside to the right of the keys of sortedData.
	evenBiggerKey = maxKey+"Z";

	ASSERT_TRUE( biggerKey.compare(evenBiggerKey) < 0);
	ASSERT_TRUE( maxKey.compare(biggerKey) < 0);
	ASSERT_TRUE( minKey.compare(maxKey) < 0);
	ASSERT_TRUE( smallerKey.compare(minKey) < 0);
	ASSERT_TRUE( evenSmallerKey.compare(smallerKey) < 0);

	// ask for keys outside, to the left, of keys(sortedData)
	qr.find(sortedData, evenSmallerKey, smallerKey);
	ASSERT_FALSE(qr.read(key,value));

	// ask for keys outside, to the right, of keys(sortedData)
	qr.find(sortedData, biggerKey, evenBiggerKey);
	ASSERT_FALSE(qr.read(key,value));

	// ask for reverse ordered keys
	qr.find(sortedData, maxKey, minKey);
	ASSERT_FALSE(qr.read(key,value));

	// ask for all the keys
	size_t count = 0;
	qr.find(sortedData, minKey, maxKey);
	for( ; qr.read(key,value); count++);
	ASSERT_EQ(count, totalCount);

	// find the smallest key  and the biggest, that are outside the range, so we encompass the entire sortedData.
	count = 0;
	qr.find(sortedData, evenSmallerKey, evenBiggerKey);
	for( ; qr.read(key,value); count++);
	ASSERT_EQ(count, totalCount);

	// pick two values inside the range of generated keys.
	StringMap::const_iterator it = verify.begin();
	std::string fk=it++->first,lk=it->first, swp;
	if( fk > lk) { swp = fk; fk = lk; lk = swp; } //make sure they're in order.
	qr.find(sortedData, fk, lk);
	//read sortedData([fk lk])
	while(qr.read(key,value))
	{
		ASSERT_EQ(value, verify[key]);
	}

	// make sure it works for a single key
	it = verify.begin();
	fk=it->first,lk=fk;
	qr.find(sortedData, fk, lk);
	qr.read(key, value);
	ASSERT_EQ( value, verify[fk] );
	
	// ask for keys (k1,k2) that are within the range defined by two consecutive keys (C1,C2) in the sorted set.
	// By definition, those keys (k1,k2) do not exist in the sorted set. The query should return false.
	// The trick is to find a couple of consecutive keys are distant enough from each other that we could
	// fit two other keys in between. Since we use keys of constant length, that means that we must find
	// (C1,C2) such that C1+=2 < C2. 
	bool foundConsecutiveKeys = false;
	std::string C1,C2,tempC1,k1,k2;
	qr.find(sortedData, minKey, maxKey);
	qr.read(C1, value);
	while( qr.read(C2, value) && !foundConsecutiveKeys)
	{
		k2 = inc(inc(C1));
		if( k2 < C2)
		{
			foundConsecutiveKeys = true;
			k1 = inc(C1);
		}
		else
		{
			C1 = C2;
		}
	}
	if(foundConsecutiveKeys)
	{
		qr.find(sortedData, k1, k2);
		ASSERT_FALSE(qr.read(key,value));
	}
	
}

TEST(find,inc)
{
	ASSERT_EQ( std::string("a"), inc("z"));
	ASSERT_EQ( std::string("ab"), inc("aa"));
	ASSERT_EQ( std::string("baa"), inc("azz"));
}
