#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "base/base.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"
#include "modules/io/spiral_file.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/io/spiral_file_mmap.h"

using namespace testing;

namespace {

struct my_serializable {
  my_serializable(const std::string contents, bool has_subpart = true)
      : m_contents(owned_membuf::from_str(contents, "spiral_file_test")) {
    if (has_subpart) {
      m_subpart.reset(new my_serializable("Uninitialized subpart", false /* no more subparts */));
    }
  }

  void create_spiral_file_part(const spiral_file_create_state &state) {
    state.set_version("my_serializable", my_version);
    state.create_membuf("contents", m_contents);
    m_mutable_contents = state.create_membuf("mutable", 100);

    if (m_subpart) {
      m_subpart->create_spiral_file_part(state.create_subpart("subpart"));
    }
  }
  void open_spiral_file_part(const spiral_file_open_state &state) {
    state.enforce_max_version("my_serializable", my_version);
    m_contents = state.open_membuf("contents");
    m_mutable_contents = state.open_mutable_membuf("mutable");

    if (m_subpart) {
      m_subpart->open_spiral_file_part(state.open_subpart("subpart"));
    }
  }

  static const product_version my_version;

  membuf m_contents;
  mutable_membuf m_mutable_contents;
  std::unique_ptr<my_serializable> m_subpart;
};

const product_version my_serializable::my_version{"1.2.3"};

}  // namespace

enum spiral_file_test_type { MEM_TEST, MMAP_TEST };

std::ostream &operator<<(std::ostream &os,
                         const std::tuple<spiral_file_test_type, spiral_file_options> &param) {
  os << "(";
  switch (std::get<0>(param)) {
    case MEM_TEST:
      os << "MEM_TEST";
    case MMAP_TEST:
      os << "MMAP_TEST";
  }
  os << ", " << std::get<1>(param) << "\n";
  return os;
}

class spiral_file_test
    : public testing::TestWithParam<std::tuple<spiral_file_test_type, spiral_file_options>> {
 public:
  std::string m_mmap_file = CONF_S(temp_root) + "/spiral_file_test";

  void SetUp() override { m_options = std::get<1>(GetParam()); }

  void create(my_serializable &my_part) {
    switch (std::get<0>(GetParam())) {
      case MEM_TEST: {
        auto mem = new spiral_file_create_mem();
        m_create_file.reset(mem);
        my_part.create_spiral_file_part(mem->create());
        break;
      }
      case MMAP_TEST: {
        unlink(m_mmap_file.c_str());
        auto mm = new spiral_file_create_mmap(m_mmap_file, m_options);
        m_create_file.reset(mm);
        my_part.create_spiral_file_part(mm->create());
        break;
      }
    }
  }

  void close() {
    if (m_create_file) {
      spiral_file_create_mem *mem = dynamic_cast<spiral_file_create_mem *>(m_create_file.get());
      if (mem) {
        m_encoded = mem->close();
      }
    }
    m_open_file.reset();
    m_create_file.reset();
  }

  void open(my_serializable &my_part) {
    switch (std::get<0>(GetParam())) {
      case MEM_TEST: {
        auto mem = new spiral_file_open_mem(m_encoded);
        m_open_file.reset(mem);
        my_part.open_spiral_file_part(mem->open());
        break;
      }
      case MMAP_TEST: {
        auto mm = new spiral_file_open_mmap(m_mmap_file, mmap_buffer::mode::read_write, m_options);
        m_open_file.reset(mm);
        my_part.open_spiral_file_part(mm->open());
        break;
      }
    }
  }

  spiral_file_mem_storage m_encoded;
  spiral_file_options m_options;
  std::unique_ptr<spiral_file_create> m_create_file;
  std::unique_ptr<spiral_file_open> m_open_file;
};

TEST_P(spiral_file_test, simple_archive) {
  auto check_metadata = [this]() {
    EXPECT_THAT(m_open_file->contents(),
                UnorderedElementsAre(
                    // For the file:
                    "file_info.json",
                    // For the top level part:
                    "part_info.json", "contents", "mutable",
                    // For the subpart:
                    "subpart/part_info.json", "subpart/contents", "subpart/mutable"));
    // UUID should be generated.
    EXPECT_NE(m_open_file->uuid(), "");
    EXPECT_THAT(m_open_file->file_info().command_line[0], HasSubstr("spiral_file_test"));
  };

  {
    my_serializable orig("Test contents");
    orig.m_subpart->m_contents = owned_membuf::from_str("Subpart contents", "spiral_file_test");
    create(orig);
    ASSERT_EQ(orig.m_mutable_contents.size(), 100);
    EXPECT_EQ("Test contents", orig.m_contents.str());
    memcpy(orig.m_mutable_contents.mutable_data(), "mutated", sizeof("mutated"));

    ASSERT_EQ(orig.m_subpart->m_mutable_contents.size(), 100);
    memcpy(orig.m_subpart->m_mutable_contents.mutable_data(), "mutated subpart",
           sizeof("mutated subpart"));
  }
  close();

  {
    my_serializable decoded("Wrong contents");
    open(decoded);
    EXPECT_EQ("mutated", std::string(decoded.m_mutable_contents.data()));
    EXPECT_EQ("Test contents", decoded.m_contents.str());
    memcpy(decoded.m_mutable_contents.mutable_data(), "Mutated again", sizeof("Mutated again"));
    EXPECT_EQ("mutated subpart", std::string(decoded.m_subpart->m_mutable_contents.data()));
    memcpy(decoded.m_subpart->m_mutable_contents.mutable_data(), "Mutated subpart again",
           sizeof("Mutated subpart again"));

    check_metadata();

    close();
  }

  {
    my_serializable decoded("Wrong contents");
    open(decoded);
    EXPECT_EQ("Mutated again", std::string(decoded.m_mutable_contents.data()));
    EXPECT_EQ("Mutated subpart again", std::string(decoded.m_subpart->m_mutable_contents.data()));
    EXPECT_EQ("Test contents", decoded.m_contents.str());

    check_metadata();

    close();
  }
}

std::vector<spiral_file_options> all_options() {
  spiral_file_options defaults;
  defaults.small_object_threshold = 1;
  spiral_file_options read_into_ram;
  read_into_ram.read_into_ram = true;
  read_into_ram.small_object_threshold = 1;
  spiral_file_options delayed_write;
  delayed_write.delayed_write = true;
  return {defaults, read_into_ram, delayed_write};
}

INSTANTIATE_TEST_CASE_P(mem_spiral_file_tests_tests, spiral_file_test,
                        ::testing::Combine(::testing::Values(MEM_TEST),
                                           ::testing::ValuesIn(all_options())));
INSTANTIATE_TEST_CASE_P(mmap_spiral_file_tests, spiral_file_test,
                        ::testing::Combine(::testing::Values(MMAP_TEST),
                                           ::testing::ValuesIn(all_options())));

// TODO(nils): Test error conditions better, for instance:
// * Missing version info when writing
// * Missing version check when reading
// * Version too new
// * Writing to read-only mmaps
// These probably require death tests.
