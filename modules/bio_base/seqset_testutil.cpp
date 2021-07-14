#include "modules/bio_base/seqset_testutil.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_mapred/make_readmap.h"
#include "modules/build_seqset/builder.h"
#include "modules/build_seqset/expand.h"
#include "modules/build_seqset/part_repo.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/output_stream.h"
#include "modules/test/test_utils.h"

#include <gtest/gtest.h>

unsigned g_seqset_build_partition_depth = 1;

std::unique_ptr<seqset_file> seqset_for_reads(
    const std::vector<dna_sequence>& reads) {
  static std::atomic<size_t> g_test_runs{0};
  size_t test_run = g_test_runs.fetch_add(1);

  build_seqset::part_repo entries(
      g_seqset_build_partition_depth,
      make_path("build_seqset_ref" + std::to_string(test_run)),
      make_path("build_seqset_repo" + std::to_string(test_run)));
  entries.open_write_pass("initial");
  for (const auto& r : reads) {
    entries.write(r, 1 /* fwd entries */, 1 /* reverse entries */);
  }
  entries.flush();

  {
    build_seqset::expander expand(entries, true /* keep temporary files */);
    expand.sort_and_dedup("", "initial", "init_sorted", "", 0, 0);
    expand.expand("init_sorted", "init_expanded", 16, 255);
    expand.sort_and_dedup("init_sorted", "init_expanded", "pass2_sorted",
                          "pass2_expanded", 1, 15);
    expand.sort_and_dedup("pass2_sorted", "pass2_expanded", "complete", "", 0,
                          0);
    size_t more_expand_needed =
        expand.expand("complete", "complete_expanded", 1, 255);
    size_t dedupped = expand.sort_and_dedup("complete", "complete_expanded",
                                            "test_out", "", 0, 0);
    CHECK_EQ(more_expand_needed, dedupped);
  }

  std::string seqset_path = make_path("test_seqset" + std::to_string(test_run));

  {
    spiral_file_create_mmap c(seqset_path);
    build_seqset::builder b;
    b.build_chunks(entries, "complete");
    b.make_seqset(c.create());
  }

  return make_unique<seqset_file>(seqset_path);
}

std::unique_ptr<seqset_flat> seqset_flat_for_seqset(const seqset* the_seqset) {
  spiral_file_mem_storage encoded;
  {
    spiral_file_create_mem c;
    seqset_flat_builder b(the_seqset);
    b.build(c.create());
    encoded = c.close();
  }

  spiral_file_open_mem o(encoded);
  return make_unique<seqset_flat>(o.open(), the_seqset);
}

std::unique_ptr<readmap> readmap_for_reads(
    const std::shared_ptr<seqset>& the_seqset,
    const std::vector<std::pair<dna_sequence, dna_sequence>>& paired_reads,
    const std::vector<dna_sequence>& unpaired_reads,
    boost::optional<std::string*> readmap_filename) {
  static std::atomic<size_t> g_test_runs{0};
  size_t test_run = g_test_runs.fetch_add(1);
  unsigned max_read_len = 0;

  manifest reads_manifest;
  output_stream_params osp;
  osp.encoding = "null";
  auto sink =
      osp.build(make_path("readmap_for_reads" + std::to_string(test_run)),
                "corrected_reads", reads_manifest);

  for (const auto& read : paired_reads) {
    corrected_reads cr;
    corrected_read cr1;
    cr1.corrected = cr1.sequence = read.first.as_string();
    cr.push_back(cr1);
    max_read_len = std::max<unsigned>(max_read_len, read.first.size());

    corrected_read cr2;
    cr2.corrected = read.second.as_string();
    cr.push_back(cr2);
    max_read_len = std::max<unsigned>(max_read_len, read.second.size());
    sink->write_msgpack<std::string, corrected_reads>("", cr);
  }

  for (const auto& read : unpaired_reads) {
    corrected_reads cr;
    corrected_read cr1;
    cr1.corrected = cr1.sequence = read.as_string();
    cr.push_back(cr1);
    max_read_len = std::max<unsigned>(max_read_len, read.size());
    sink->write_msgpack<std::string, corrected_reads>("", cr);
  }
  sink->close();
  sink.reset();

  std::string readmap_file_path =
      make_path("test_readmap" + std::to_string(test_run));
  if (readmap_filename) {
    **readmap_filename = readmap_file_path;
  }
  make_readmap::do_make(readmap_file_path, *the_seqset, reads_manifest,
                        true /* is paired */, max_read_len);

  return make_unique<readmap>(the_seqset, readmap_file_path);
}

void dump_seqset(const std::string& prefix, const seqset& the_seqset) {
  for (size_t i = 0; i < the_seqset.size(); ++i) {
    std::cerr << prefix << i << ": "
              << the_seqset.ctx_entry(i).sequence().as_string() << "\n";
  }
}

std::pair<std::shared_ptr<seqset_file>, std::unique_ptr<readmap>>
biograph_for_reads(const std::vector<std::vector<dna_sequence>>& all_reads) {
  std::vector<dna_sequence> seq_input;
  std::vector<std::pair<dna_sequence, dna_sequence>> paired;
  std::vector<dna_sequence> unpaired;
  for (auto i : all_reads) {
    seq_input.push_back(i[0]);
    if (i.size() == 2 && i[1].size() != 0) {
      seq_input.push_back(i[1]);
      paired.push_back(std::make_pair(i[0], i[1]));
    } else {
      unpaired.push_back(i[0]);
    }
  }
  std::shared_ptr<seqset> new_seqset = seqset_for_reads(seq_input);
  auto new_readmap = readmap_for_reads(new_seqset, paired, unpaired);
  return std::make_pair< std::shared_ptr<seqset>,
                         std::unique_ptr<readmap> >(std::move(new_seqset),
                                  std::move(new_readmap));
}

void PrintTo(const seqset_range& r, std::ostream* os) {
  if (r.valid()) {
    (*os) << "[" << r.begin() << "-" << r.end() << ") " << r.sequence();
  } else {
    (*os) << "[invalid seqset_range]";
  }
}
