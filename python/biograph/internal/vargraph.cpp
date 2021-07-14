#include "python/biograph/internal/vargraph.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"

#include <boost/algorithm/string.hpp>


py_list coverage_wrapper::get_base_cov() {
  py_list ret;
  for(auto iter = m_cov.base_cov.begin(); iter != m_cov.base_cov.end(); iter++)
  {
    ret.append(*iter);
  }
  return ret;
}

py_list coverage_wrapper::get_span_cov() {
  py_list ret;
  for(auto iter = m_cov.span_cov.begin(); iter != m_cov.span_cov.end(); iter++)
  {
    ret.append(*iter);
  }
  return ret;
}

py_list node_wrapper::get_upstream_edges() {
  py_list ret;
  for (const vargraph::edge_t* e : m_node.upstream) {
    ret.append(edge_wrapper(m_graph, *e));
  }
  return ret;
}

py_list node_wrapper::get_downstream_edges() {
  py_list ret;
  for (const vargraph::edge_t* e : m_node.downstream) {
    ret.append(edge_wrapper(m_graph, *e));
  }
  return ret;
}

vargraph_wrapper::vargraph_wrapper(const dna_sequence& contig, size_t min_pair, size_t max_pair)
    : m_graph(std::make_shared<vargraph>(contig, min_pair, max_pair)) 
    , m_contig_size(contig.size())
    , m_start(0)
    , m_end(contig.size())
{}

//I think this not being a slice could be the problem?
vargraph_wrapper::vargraph_wrapper(const reference_range& ref, size_t min_pair, size_t max_pair)
    : m_graph(std::make_shared<vargraph>(ref.sequence(), min_pair, max_pair)) 
    , m_contig_size(ref.size())
    , m_start(ref.get_start())
    , m_end(ref.get_end())
{}

void vargraph_wrapper::add_variant(uint32_t start, uint32_t end, const dna_sequence& seq) {
  m_graph->add_variant(start, end, seq);
}

void vargraph_wrapper::trace(const std::shared_ptr<biograph> bg) {
  trace_sub(bg, 0, int(m_contig_size));
}

void vargraph_wrapper::trace_sub(const std::shared_ptr<biograph> bg, uint32_t start, uint32_t end) {
  auto ss = bg->get_seqset();   //->get_seqset();
  auto rm = bg->open_readmap("");  // or whichever one is active
  m_graph->trace(ss->get_seqset(), *rm, start, end);
} 

/* void vargraph_wrapper::trace_sub(const std::shared_ptr<biograph> bg, uint32_t start, uint32_t end) {
  printf("Hello world\n");
  Config::load("/home/english/spiral/etc/products/unittest.json");
  printf("Loading reference\n");
  Config::set("reference_path", "/share/reference/hs37d5");
  reference ref("");
  printf("Loading seqset\n");
  seqset_file file("/mnt/data/ajtrio/manual_builder/37/results/HG002.bg/seqset");
  const seqset& ss = file.get_seqset();
  printf("Loading readmap\n");
  readmap rm(ss, "/mnt/data/ajtrio/manual_builder/37/results/HG002.bg/coverage/cf3236d07a8a8d22a7784274c2173036939515c1.readmap");
  // Single supercontig
  std::string chr = "1";
  start = 57238525 + 943;
  end = 57252700;
  size_t flat_start = ref.get_assembly().flatten(chr, start);
  size_t flat_end = ref.get_assembly().flatten(chr, end);
  dna_slice slice(ref.get_dna(flat_start), ref.get_dna(flat_end));
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
  vg->trace(ss, rm, 0, slice.size());
  printf("Done, time = %ld\n", time(0));
  std::string ret = ""; //Dump of vargraph with full_cov = " + full_cov ? "true" : "false";
  for(const auto& kvp : vg->get_nodes()) {
    ret += kvp.second->as_string() + "\n";
  }
  for(const auto& e : vg->get_edges()) {
    ret += e->upstream->as_string() + " ->" + e->downstream->as_string() + "\n";
    ret += "  unpaired: " + std::to_string(e->unpaired) + ", paired: " + std::to_string(e->paired) + "\n";
  }
  std::cout << ret;
} */

// should have a way to get a single position
// py_list vargraph_wrapper::get_nodes(uint32_t pos)
py_list vargraph_wrapper::get_nodes() {
  py_list ret;
  for (const auto& kvp: m_graph->get_nodes()) {
    ret.append(node_wrapper(m_graph, *kvp.second));
  }
  return ret;
}

py_list vargraph_wrapper::get_edges() {
  py_list ret;
  for (const vargraph::edge_t* e : m_graph->get_edges()) {
    ret.append(edge_wrapper(m_graph, *e));
  }
  return ret;
}

std::string vargraph_wrapper::dump_graph() {
  std::string ret = ""; //Dump of vargraph with full_cov = " + full_cov ? "true" : "false";
  for(const auto& kvp : m_graph->get_nodes()) {
    ret += kvp.second->as_string() + "\n";
  }
  for(const auto& e : m_graph->get_edges()) {
    ret += e->upstream->as_string() + " ->" + e->downstream->as_string() + "\n";
    ret += "  unpaired: " + std::to_string(e->unpaired) + ", paired: " + std::to_string(e->paired) + "\n";
  }
  return ret;
}

using namespace pybind11;

void bind_vargraph(module& m) {
  class_<vargraph_wrapper>(m, "VarGraph", "Basic Block\n")
      .def(init<const dna_sequence&>())
      .def(init<const dna_sequence&, size_t, size_t>())
      .def("add_variant", &vargraph_wrapper::add_variant, "Doc\n")
      .def("trace", &vargraph_wrapper::trace, "Doc\n")
      .def("trace_sub", &vargraph_wrapper::trace_sub, "Doc\n")
      .def("get_nodes", &vargraph_wrapper::get_nodes, "Doc\n")
      .def("get_edges", &vargraph_wrapper::get_edges, "Doc\n")
      .def("dump_graph", &vargraph_wrapper::dump_graph, "Outputs all nodes/edges\n");

  class_<node_wrapper>(m,
    "VarNode", "Genomic sequence\n")
    .def_property_readonly("is_ref", &node_wrapper::get_is_ref)
    .def_property_readonly("start", &node_wrapper::get_start)
    .def_property_readonly("end", &node_wrapper::get_end)
    .def_property_readonly("seq", &node_wrapper::get_seq)
    .def_property_readonly("paired", &node_wrapper::get_paired)
    .def_property_readonly("unpaired", &node_wrapper::get_unpaired)
    .def_property_readonly("upstream", &node_wrapper::get_upstream_edges)
    .def_property_readonly("downstream", &node_wrapper::get_downstream_edges)
    .def("__repr__", &node_wrapper::repr)
  ;

  class_<coverage_wrapper>(m,
    "VarCoverage", "Coverage information for a VarNode\n")
    .def_property_readonly("base_cov", &coverage_wrapper::get_base_cov,
         "Doc \n")
    .def_property_readonly("span_cov", &coverage_wrapper::get_span_cov,
         "Doc\n")
  ;

  class_<edge_wrapper>(m,
    "VarEdge", "Edge between two VarNodes\n")
    .def_property_readonly("up_node", &edge_wrapper::get_upstream_edge,
        "Doc\n")
    .def_property_readonly("dn_node", &edge_wrapper::get_downstream_edge,
        "Doc\n")
    .def_property_readonly("paired", &edge_wrapper::get_paired,
        "Doc\n")
    .def_property_readonly("unpaired", &edge_wrapper::get_unpaired,
        "Doc\n")
    .def("__repr__", &edge_wrapper::repr)
  ;
}
