#include "modules/bio_base/make_mergemap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/bio_base/seqset_merger.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file_mmap.h"

#include <boost/filesystem.hpp>
#include <vector>

static void update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

int main(int argc, char** argv) {
  spiral_init(&argc, &argv);
  log_init(nullptr, 2, true);

  SPLOG("Logging test");

  std::vector<std::string> inputs = {
      "SRR2098210.seqset", "SRR2098211.seqset"
      //"SRR2098212.seqset",
      //"SRR2098213.seqset",
  };

  for (auto& input : inputs) {
    input = "/scratch/" + input;
  }

  for (const auto& input : inputs) {
    std::string path = input + ".flat";
    if (boost::filesystem::exists(path)) {
      SPLOG("Skipping flat generation; %s already exists", path.c_str());
      continue;
    }

    SPLOG("Building flat seqset for %s", input.c_str());
    std::unique_ptr<seqset_file> ss_f(new seqset_file(input));
    SPLOG("Creating spiral file");
    std::unique_ptr<spiral_file_create_mmap> c(
        new spiral_file_create_mmap(path));
    SPLOG("Making pop front cache");
    ss_f->get_seqset().populate_pop_front_cache(update_progress);
    SPLOG("Creating flat output");
    std::unique_ptr<seqset_flat_builder> flat(
        new seqset_flat_builder(&ss_f->get_seqset()));
    SPLOG("Building flat");
    flat->build(c->create(), update_progress);
    SPLOG("Flat build complete");
    flat.reset();
    c->close();
    ss_f->get_seqset().clear_pop_front_cache();
  }

  SPLOG("Opening flats");
  std::vector<std::unique_ptr<seqset_file>> seqsets;
  std::vector<std::unique_ptr<seqset_flat>> flats;
  std::vector<const seqset_flat*> flat_ptrs;

  for (const auto& input : inputs) {
    seqsets.emplace_back(new seqset_file(input));
    spiral_file_open_mmap o(input + ".flat");
    flats.emplace_back(
        new seqset_flat(o.open(), &seqsets.back()->get_seqset()));
    flat_ptrs.emplace_back(flats.back().get());
  }

  bool all_mergemaps_done = true;
  for (const auto& input : inputs) {
    std::string path = input + ".mergemap";
    if (!boost::filesystem::exists(path)) {
      all_mergemaps_done = false;
    }
  }

  unlink("/scratch/merged.seqset");
  spiral_file_create_mmap create_merge("/scratch/merged.seqset");

  if (all_mergemaps_done) {
    SPLOG("All mergemaps already done; skipping mergemap generation");
  } else {
    SPLOG("Building mergemaps");
    make_mergemap mm_make(flat_ptrs);
    mm_make.build(update_progress);

    SPLOG("%lu entries in resultant merge; writing mergemaps",
          mm_make.total_merged_entries());
    unsigned input_index = 0;
    for (const auto& input : inputs) {
      std::string path = input + ".mergemap";
      unlink(path.c_str());
      spiral_file_create_mmap c(path);
      seqset_mergemap_builder build_mergemap(
          c.create(), seqsets[input_index]->get_seqset().uuid(),
          create_merge.uuid(), mm_make.total_merged_entries());
      mm_make.fill_mergemap(input_index++, &build_mergemap, update_progress);
      c.close();
    }
  }

  SPLOG("Opening mergemaps");
  std::vector<std::unique_ptr<seqset_mergemap>> mergemaps;
  std::vector<const seqset_mergemap*> mergemap_ptrs;
  for (const auto& input : inputs) {
    spiral_file_open_mmap o(input + ".mergemap");
    mergemaps.emplace_back(new seqset_mergemap(o.open()));
    mergemap_ptrs.emplace_back(mergemaps.back().get());
  }

  SPLOG("Generting final merge.");
  seqset_merger merger(flat_ptrs, mergemap_ptrs);
  merger.build(create_merge.create(), update_progress);
  create_merge.close();
  SPLOG("All done");
}
