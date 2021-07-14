#pragma once

#include <pybind11/pybind11.h>

#include "modules/vargraph/vargraph.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

typedef pybind11::list py_list;
typedef pybind11::dict py_dict;

class edge_wrapper;

class coverage_wrapper {
 public:
  coverage_wrapper(const std::shared_ptr<vargraph>& graph, const vargraph::cov_info_t& cov)
      : m_graph(graph), m_cov(cov) {}
  py_list get_base_cov();
  py_list get_span_cov();
 private:
  std::shared_ptr<vargraph> m_graph;
  const vargraph::cov_info_t& m_cov;
};

class node_wrapper {
 public:
  node_wrapper(const std::shared_ptr<vargraph>& graph, const vargraph::node_t& node)
      : m_graph(graph), m_node(node) {}

  bool get_is_ref() { return m_node.is_ref; }
  uint32_t get_start() { return m_node.start; }
  uint32_t get_end() { return m_node.end; }
  dna_sequence get_seq() { return dna_sequence(m_node.seq); }
  coverage_wrapper get_paired() { return coverage_wrapper(m_graph, m_node.paired); }
  coverage_wrapper get_unpaired() { return coverage_wrapper(m_graph, m_node.unpaired); }

  py_list get_upstream_edges();
  py_list get_downstream_edges();
  std::string repr() { return "<VarNode [" + std::to_string(m_node.start) + ", " + std::to_string(m_node.end) + 
                              ") is_ref=" + std::to_string(m_node.is_ref) + " sz=" + std::to_string(m_node.seq.size()) 
                              + ">"; }
 private:
  std::shared_ptr<vargraph> m_graph;
  const vargraph::node_t& m_node;
};

class edge_wrapper {
 public:
  edge_wrapper(const std::shared_ptr<vargraph>& graph, const vargraph::edge_t& edge)
      : m_graph(graph), m_edge(edge) {}

  node_wrapper get_upstream_edge() { return node_wrapper(m_graph, *m_edge.upstream); }
  node_wrapper get_downstream_edge() { return node_wrapper(m_graph, *m_edge.downstream); }

  uint32_t get_paired() { return m_edge.paired; }
  uint32_t get_unpaired() { return m_edge.unpaired; }

  std::string repr() { return "<VarEdge (" + std::to_string(get_upstream_edge().get_end()) + ") up_is_ref=" + 
                              std::string(get_upstream_edge().get_is_ref() ? "true" : "false") + " dn_is_ref=" + 
                              std::string(get_downstream_edge().get_is_ref() ? "true" : "false") + ">"; }
 private:
  std::shared_ptr<vargraph> m_graph;
  const vargraph::edge_t& m_edge;
};

class vargraph_wrapper {
 public:
  // Constructor -- can't use dna_slice 
  vargraph_wrapper(const dna_sequence& contig, size_t min_pair = 100, size_t max_pair = 1000);
  // Constructor -- reference range
  vargraph_wrapper(const reference_range& ref, size_t min_pair = 100, size_t max_pair = 1000);

  // Accessors 
  void add_variant(uint32_t start, uint32_t end, const dna_sequence& seq);

  // How do I call this multiple times so I can trace PER readmap?
  void trace(const std::shared_ptr<biograph> bg);
  void trace_sub(const std::shared_ptr<biograph> bg, uint32_t start, uint32_t end);

  //This will take some work to collapse a location
  // to allow per-allele counting
  // returns a dict - {pos:[node, ..], ..}
  py_list get_nodes();
  // returns a list - [(edge, edge), ..}
  py_list get_edges();

  std::string dump_graph();
 private:
  std::shared_ptr<vargraph> m_graph;
  size_t m_contig_size;
  uint32_t m_start; // absolute start position of the sequence
  uint32_t m_end;   // absolute end position of the sequence

};

