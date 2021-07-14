
#include "vargraph.h"

#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"

#include <boost/algorithm/string.hpp>

std::string str_vec(const std::vector<int>& vec) {
  std::string r = "[";
  for(size_t i = 0; i < vec.size(); i++) {
    r += std::to_string(vec[i]);
    if (i != vec.size() - 1) { r += " "; }
  }
  r += "]";
  return r;
}

void dump_cov(const char* name, const vargraph::cov_info_t& ci) {
  printf("  %s:\n", name);
  printf("    base_cov: %s\n", str_vec(ci.base_cov).c_str());
  printf("    span_cov: %s\n", str_vec(ci.span_cov).c_str());
};

void dump_graph(vargraph* vg, bool full_cov = true)
{
  printf("Dump of vargraph with full_cov = %s\n", full_cov ? "true" : "false");
  for(const auto& kvp : vg->get_nodes()) {
    printf("%s\n", kvp.second->as_string().c_str());
    if (full_cov) {
      dump_cov("unpaired", kvp.second->unpaired);
      dump_cov("paired", kvp.second->paired);
    }
  }
  for(const auto& e : vg->get_edges()) {
    printf("%s->%s\n", e->upstream->as_string().c_str(), e->downstream->as_string().c_str());
    printf("  unpaired: %u, paired: %u\n", e->unpaired, e->paired);
  }
}

void dump_reads(std::vector<std::vector<dna_sequence>> fake_reads) {
  for(auto const& read : fake_reads) {
    if (read.size() == 1) {
      printf("unpr %s\n", read[0].as_string().c_str());
    } else if (read.size() == 2) {
      printf("pair %s %s\n", read[0].as_string().c_str(), read[1].as_string().c_str());
    }
  }
}


int main() {
  printf("Hello world\n");
  Config::load("/home/english/spiral/etc/products/unittest.json");
  printf("Loading reference\n");
  Config::set("reference_path", "/share/reference/hs37d5");
  reference ref("");
  printf("Loading seqset\n");
  auto ss = std::make_shared<seqset>("/mnt/data/ajtrio/manual_builder/37/results/HG002.bg/seqset");
  printf("Loading readmap\n");
  readmap rm(ss, "/mnt/data/ajtrio/manual_builder/37/results/HG002.bg/coverage/cf3236d07a8a8d22a7784274c2173036939515c1.readmap");
  // Single supercontig
  std::string chr = "1";
  size_t start = 37313385;
  size_t end = 37322046;
  size_t flat_start = ref.get_assembly().flatten(chr, start);
  size_t flat_end = ref.get_assembly().flatten(chr, end);
  //dna_slice slice(ref.get_dna(flat_start), ref.get_dna(flat_end));
  dna_sequence slice(ref.get_dna(flat_start), ref.get_dna(flat_end));
  printf("Got a slice, size = %d\n", int(slice.size()));
  auto vg = std::make_shared<vargraph>(slice);
  file_reader fr("/home/english/data/single_sample_pcmp/debugging/1:37248525-77248525_calls.vcf");
  printf("Loading VCF, time = %ld\n", time(0));
  std::string line;
  int cnt = 0;
  while (fr.readline(line, 500000)) {
    std::vector<std::string> fields;
    if (line[0] == '#') continue;
    cnt += 1;
    boost::split(fields, line, boost::is_any_of("\t"));
    long pos = atol(fields[1].c_str());
    std::vector<std::string> alts;
    boost::split(alts, fields[4], boost::is_any_of(","));
    for(const auto& alt : alts) {
      vg->add_variant(pos - start - 1, pos - start - 1 + fields[3].size(), dna_sequence(alt));
      printf("Adding variant, [%lu-%lu) local coords, alt = %s\n", pos - start - 1, pos - start - 1 + fields[3].size(), alt.c_str());
    }
  }
  printf("Found %d\n", cnt);
  printf("Doing trace, %d %d @ time = %ld\n", 0, int(slice.size()), time(0));
  vg->trace(*ss, rm, 0, slice.size());
  printf("Done, time = %ld\n", time(0));
  dump_graph(vg.get(), false);
}
