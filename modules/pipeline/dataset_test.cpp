#include "modules/pipeline/dataset_path.h"
#include "modules/pipeline/primitives.h"
#include "modules/pipeline/dataset_meta.h"
#include "gtest/gtest.h"
#include "modules/web/couchdb.h"
#include "modules/io/config.h" 
#include "modules/pipeline/ottoman.h"

class dataset_test : public ::testing::Test {
public:
	void SetUp() override {
		if (!os) {
			os.reset(new ottoman_server());
		}
	}
private:
	// Dataset tests depend on ottoman being accessible, but ottoman
	// server doesn't like being taken down.
	static std::unique_ptr<ottoman_server> os;
};
std::unique_ptr<ottoman_server> dataset_test::os;

TEST_F(dataset_test, basic)
{
	dataset_path p("/api/users/spiral_tester/data/test_a");

	ASSERT_EQ(p.friendly().compare("/test_a"), 0);
	ASSERT_EQ(p.url().compare("/api/users/spiral_tester/data/test_a"), 0);
	ASSERT_EQ(p.user().compare("spiral_tester"), 0);
	ASSERT_EQ(p.parent().compare("/api/users/spiral_tester/data"), 0);
	ASSERT_EQ(p.name().compare("test_a"), 0);
	ASSERT_EQ(p.root().append(p.name()).url(), p.url());
	ASSERT_FALSE(p.is_reference());
}

void check_ls(const std::string& path, size_t expected_size)
{
	dataset_path dp(path);
	std::vector<direntry> listing = dp.list_dir();
	ASSERT_EQ(listing.size(), expected_size);
}

void recursive_rmdir(const dataset_path& dp) 
{
	auto exists = dp.exists();
	if (exists == path::e_directory) {
		std::vector<direntry> listing = dp.list_dir();
		for (const auto& item : listing) {
			dataset_path child(item.url);
			recursive_rmdir(child);
		}
		dp.remove();
	} else if (exists == path::e_file) {
		dp.remove();
	}
}

TEST_F(dataset_test, need_couchdb)
{
	gen_cache(nullptr);

	dataset_path root("/api/users/spiral_tester/data");
	recursive_rmdir(root);

	dataset_path p("/api/users/spiral_tester/data/test_a/my_dir/mo_dir");

	ASSERT_NO_THROW(p.remove());
	ASSERT_EQ(p.exists(), path::e_no_exist);
	check_ls("/api/users/spiral_tester/data", 0);
	check_ls("/api/users/spiral_tester/data/test_a", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir/mo_dir", 0);

	// triggers the recursive mkdirs for test_a , my_dir and mo_dir
	ASSERT_NO_THROW(p.mkdir()); 

	check_ls("/api/users/spiral_tester/data", 1);
	check_ls("/api/users/spiral_tester/data/test_a", 1);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir", 1);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir/mo_dir", 0);

	dataset_path p_dir("/api/users/spiral_tester/data/test_a");
	dataset_path p_my_dir(p_dir.append("my_dir"));

	ASSERT_EQ(p_dir.exists(), path::e_directory);
	ASSERT_EQ(p_my_dir.exists(), path::e_directory);

	ASSERT_NO_THROW(p.remove());
	check_ls("/api/users/spiral_tester/data", 1);
	check_ls("/api/users/spiral_tester/data/test_a", 1);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir/mo_dir", 0);

	ASSERT_EQ(p_my_dir.exists(), path::e_directory);

	ASSERT_NO_THROW(p_my_dir.remove());
	check_ls("/api/users/spiral_tester/data", 1);
	check_ls("/api/users/spiral_tester/data/test_a", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir/mo_dir", 0);

	ASSERT_NO_THROW(p_dir.remove());
	check_ls("/api/users/spiral_tester/data", 0);
	check_ls("/api/users/spiral_tester/data/test_a", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir", 0);
	check_ls("/api/users/spiral_tester/data/test_a/my_dir/mo_dir", 0);

	ASSERT_EQ(p_dir.exists(), path::e_no_exist);
	ASSERT_EQ(p_my_dir.exists(), path::e_no_exist);
	ASSERT_EQ(p.exists(), path::e_no_exist);

	dataset_path p2("/api/users/spiral_tester/data/my_file");
	dataset_meta dm;
	dm.type = datatype_registry::find("unaligned_reads");

	ASSERT_NO_THROW(p2.remove());
	ASSERT_NO_THROW(p2.create(dm));
	ASSERT_THROW(p2.mkdir(), io_exception);
	dataset_path p3(p2.append("dir3"));
	ASSERT_ANY_THROW(p3.mkdir());

	dataset_path p4("/api/users/spiral_tester/data");
	std::vector<direntry> listing = p4.list_dir();
	ASSERT_EQ(listing.size(), (size_t)1);
	if (listing.size() == 1) {
		ASSERT_EQ(listing[0]._id, "/api/users/spiral_tester/data/my_file");
	}

	dataset_path p5("/api/users/spiral_tester/data/my_dir");
	ASSERT_NO_THROW(p5.mkdir()); 

	listing = p4.list_dir();
	ASSERT_EQ(listing.size(), (size_t)2);
	if (listing.size() == 2) {
		ASSERT_EQ(listing[0]._id, "/api/users/spiral_tester/data/my_dir");
		ASSERT_EQ(listing[1]._id, "/api/users/spiral_tester/data/my_file");
	}

	auto foo = p5.append("foo");
	foo.mkdir();
	ASSERT_NO_THROW( p4.remove(true) );
	ASSERT_EQ( path::e_no_exist, foo.exists() );
	ASSERT_EQ( path::e_no_exist, p2.exists() );
	ASSERT_EQ( path::e_no_exist, p5.exists() );
	ASSERT_EQ( path::e_no_exist, p4.exists() );
}
