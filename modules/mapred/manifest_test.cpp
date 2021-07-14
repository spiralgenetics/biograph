#include <atomic>

#include "modules/io/utils.h"
#include "modules/io/file_io.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/manifest_parallel.h"
#include "modules/bio_base/corrected_read.h"
#include "gtest/gtest.h"

class int_sorter : public sorter
{
public:
	int_sorter(const std::string& str) 
	{}

	int compare(const std::string& key1, const std::string& key2) const override
	{
		int k1 = atoi(key1.c_str());
		int k2 = atoi(key2.c_str());
		if (k1 < k2) {
			return -1;
		}
		if (k2 < k1) {
			return 1;
		}
		return 0;
	}

	std::string bump_back(const std::string& key) const override
	{ 
		return printstring("%d", std::max(0, atoi(key.c_str()) - 100)); 
	}

	size_t partition(const std::string& key, size_t num_partitions) const override
	{ 
		return atoi(key.c_str()) % num_partitions; 
	}
};

REGISTER_1(sorter, int, const std::string&);

TEST(manifest, int_sort_reduce)
{
	manifest m("int");
	m.add(file_info(path("f1"), 1000, 10, "0", "20"), 0);
	m.add(file_info(path("f2"), 1000, 10, "20", "130"), 0);
	m.add(file_info(path("f3"), 1000, 10, "130", "150"), 0);
	m.add(file_info(path("f4"), 1000, 10, "152", "300"), 0);
	m.add(file_info(path("f5"), 1000, 10, "301", "330"), 0);
	m.add(file_info(path("f6"), 1000, 10, "330", "380"), 0);
	m.add(file_info(path("f7"), 1000, 10, "381", "430"), 0);
	m.add(file_info(path("f8"), 1000, 10, "435", "550"), 0);
	m.add(file_info(path("f9"), 1000, 10, "551", "623"), 0);
	m.add(file_info(path("f10"), 1000, 10, "631", "831"), 0);
	m.add(file_info(path("f11"), 1000, 10, "840", "950"), 0);
	m.add(file_info(path("f12"), 1000, 10, "950", "960"), 0);
	m.add(file_info(path("f13"), 1000, 10, "960", "970"), 0);
	m.add(file_info(path("f14"), 1000, 10, "972", "1020"), 0);
	m.add(file_info(path("f15"), 1000, 10, "1020", "1050"), 0);
	m.add(file_info(path("f16"), 1000, 10, "1050", "1080"), 0);
	m.add(file_info(path("f17"), 1000, 10, "1080", "1110"), 0);
	m.add(file_info(path("f18"), 1000, 10, "1110", "1140"), 0);
	m.add(file_info(path("f19"), 1000, 10, "1140", "1170"), 0);
	m.add(file_info(path("f20"), 1000, 10, "1170", "1200"), 0);
	m.add(file_info(path("f21"), 1000, 10, "1200", "1230"), 0);
	m.add(file_info(path("f22"), 1000, 10, "1230", "1260"), 0);
	m.add(file_info(path("f23"), 1000, 10, "1260", "1290"), 0);
	m.add(file_info(path("f24"), 1000, 10, "1290", "1320"), 0);
	m.add(file_info(path("f25"), 1000, 10, "1500", "1523"), 0);
	m.add(file_info(path("f26"), 1000, 10, "1523", "1700"), 0);
	m.add(file_info(path("f27"), 1000, 10, "1705", "2000"), 0);
	std::vector<input_stream_params> split;
	m.split_sort_reduce(split, 2500, false);
	for (size_t i = 0; i < split.size(); i++) {
		const input_stream_params& is = split[i];
		printf("From: '%s' -> '%s'\n", is.begin_on.c_str(), is.end_before.c_str());
		for (size_t j = 0; j < is.inputs.size(); j++) {
			const file_info& fi = is.inputs[j];
			printf("   %s: '%s'-'%s'\n", fi.file.url().c_str(), fi.first_key.c_str(), fi.last_key.c_str());
		}
	}
	ASSERT_EQ(split.size(), size_t(7));
}


TEST(manifest, split_sort_reduce)
{
	manifest m("lexical");
	m.add(file_info(path("f1"), 1000, 10, "A", "B"), 0);
	m.add(file_info(path("f2"), 1000, 10, "B", "B"), 0);
	m.add(file_info(path("f3"), 1000, 10, "B", "B"), 0);
	m.add(file_info(path("f4"), 1000, 10, "B", "C"), 0);
	m.add(file_info(path("f5"), 1000, 10, "C", "F"), 0);
	m.add(file_info(path("f6"), 1000, 10, "G", "J"), 0);
	m.add(file_info(path("f7"), 1000, 10, "J", "K"), 0);
	m.add(file_info(path("f8"), 1000, 10, "K", "N"), 0);
	m.add(file_info(path("f9"), 1000, 10, "N", "O"), 0);
	m.add(file_info(path("f10"), 1000, 10, "P", "R"), 0);
	m.add(file_info(path("f11"), 1000, 10, "S", "U"), 0);
	m.add(file_info(path("f12"), 1000, 10, "U", "X"), 0);
	m.add(file_info(path("f13"), 1000, 10, "X", "X"), 0);
	m.add(file_info(path("f14"), 1000, 10, "X", "X"), 0);
	m.add(file_info(path("f15"), 1000, 10, "X", "X"), 0);
	m.add(file_info(path("f16"), 1000, 10, "X", "X"), 0);
	m.add(file_info(path("f17"), 1000, 10, "X", "Z"), 0);
	std::vector<input_stream_params> split;
	m.split_sort_reduce(split, 2500);
	for (size_t i = 0; i < split.size(); i++) {
		const input_stream_params& is = split[i];
		printf("From: '%s' -> '%s'\n", is.begin_on.c_str(), is.end_before.c_str());
		for (size_t j = 0; j < is.inputs.size(); j++) {
			const file_info& fi = is.inputs[j];
			printf("   %s: '%s'-'%s'\n", fi.file.url().c_str(), fi.first_key.c_str(), fi.last_key.c_str());
		}
	}
	ASSERT_EQ(split.size(), size_t(5));
}

TEST(manifest, add_manifest)
{
	// merge manifests with no encoding tag
	manifest m00("lexical");
	manifest m01("lexical");
	manifest m02("lexical");

	ASSERT_NO_THROW(m00.add(m01));
	ASSERT_NO_THROW(m00.add(m02));
}

TEST(manifest, meta)
{
	manifest m;

	int a_value = 253;
	ASSERT_NO_THROW(m.metadata().set("foo", "bar", a_value));
	ASSERT_EQ(m.metadata().get<int>("foo", "bar"), a_value);
	ASSERT_EQ(m.metadata().get("foo", "yo", 11), 11);

	ASSERT_TRUE(m.metadata().has_key("foo", "bar"));
	ASSERT_THROW(m.metadata().get<int>("", "bar"), io_exception);
	ASSERT_THROW(m.metadata().get<int>("foo", ""), io_exception);
	ASSERT_EQ(m.metadata().get<int>("foo", "bad_key", 13), 13);
	ASSERT_FALSE(m.metadata().has_key("foo", "bad_key"));
	ASSERT_EQ(m.metadata().get<int>("bad_ns", "bar", 17), 17);

	manifest m1;
	ASSERT_NO_THROW(m1.merge_tags(m));
	ASSERT_EQ(m1.metadata().get<int>("foo", "bar"), a_value);
	ASSERT_NO_THROW(m1.metadata().unset("foo", "bar"));
	ASSERT_THROW(m1.metadata().get<int>("foo", "bar"), io_exception);
}

TEST(manifest, parallel)
{
	class record_counter
	{
	public:
		record_counter() : m_count(0) {}
		record_counter(record_counter&& rhs) { std::atomic_exchange(&m_count, rhs.m_count.load()); }
		void operator()(
			const std::string& /*read_id*/
			, const corrected_reads& /*the_reads*/
			, int /* file_info_id */
			, int /*cumulative_record*/
		)
		{
			m_count++;
		}
		int get_count() const { return m_count.load(); }
	
	private:
		std::atomic<int> m_count;
	};

	std::string serialized_manifest{ slurp_file("datasets/hiv/corrected/ERR381524.corrected_reads") };
	manifest corrected_read_manifest = inline_json_deserialize<manifest>(serialized_manifest);
	auto returned_functor = manifest_parallelize<record_counter, std::string, corrected_reads>(corrected_read_manifest, record_counter());
	ASSERT_EQ(returned_functor.get_count(), corrected_read_manifest.get_num_records());
}
