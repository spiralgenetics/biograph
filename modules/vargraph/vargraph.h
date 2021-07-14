#pragma once

#include <deque>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"

// A vargraph is specific to a given supercontig YO...
class vargraph {
 public:
  struct node_t;

  struct edge_t {
    edge_t(node_t* _upstream, node_t* _downstream) : upstream(_upstream), downstream(_downstream) {}
    void flatten();  // Simplify coverage numbers
    node_t* upstream;  // Node upstream of this edege
    node_t* downstream;  // Node downstream of this edge
    // For a given read id, for a given point inside the read, is it paired?
    std::unordered_map<uint32_t, std::unordered_map<unsigned, bool>> coverage;  
    int unpaired;  // Number of unpaired reads that traverse this edge
    int paired;  // Number of paired reads that traverse this edge
    std::string as_string() const;
  };

  struct cov_info_t {
    // For each base, how many reads cover that base?
    std::vector<int> base_cov;
    // For each pair of bases (numbered by the first base), how many reads span across the two bases
    // span_cov.size() == base_cov.size() - 1
    std::vector<int> span_cov;
  };

  // A node is a sequence in genome with [start, end) relative to reference
  struct node_t {
    bool is_ref;  // Is this node reference
    uint32_t start;  // Start (inclusive) relative to reference
    uint32_t end;  // End (exclusive) relative to reference
    dna_slice seq;  // The sequence (never empty) of this node
    dna_sequence backing;  // A backing sequence (for non-ref cases)
    std::vector<edge_t*> upstream; // Upstream edges
    std::vector<edge_t*> downstream; // Downstream edges
    cov_info_t unpaired;  // Coverage for unpaired reads
    cov_info_t paired;  // Coverage for paired reads
    std::string as_string() const;
  };

 private:
  // Size of entire contig
  size_t m_contig_size; 
  
  // Initial reference sequence
  dna_sequence m_contig;

  // All node_t by start point
  std::multimap<uint32_t, node_t*> m_gnodes;

  // All edges (for memory management only)
  std::set<edge_t*> m_edges;

  // Min/Max pairing distance to consider
  size_t m_min_pair;
  size_t m_max_pair;

  // Reference node_t by start point
  std::map<uint32_t, node_t*> m_ref_nodes;

  // Coverage of a node_t by a read
  // So essentially a mapping of a read to a node
  struct read_aln_t {
    // Makes a coverage entry, and chops the end from read size
    static read_aln_t make_aln(unsigned& read_size, uint32_t pos);

    // Equality and hash functions for use in unordered maps
    bool operator==(const read_aln_t& rhs) const;
    struct hash_t {
      size_t operator()(const read_aln_t& x) const;
    };

    // Position of alignment in terms of both read (query) position, as well
    // as within the node (query).  Both numbers are 0-based indexes into local
    // coordinates
    uint32_t target_start;  // The start of coverage in the entry -target start
    uint32_t target_end;  // The end of coverage in the entry -target end
    unsigned query_start;  // The start of the read's coverage -query start
    unsigned query_end;    // The end of the read's coverage -query end
  };

  // For each node, all the reads that cover it and if they are paired
  struct aln_info_t {
    void flatten(node_t* cur_node) const;
    // { read_id -> { read_aln -> is_paired } }
    std::unordered_map<uint32_t, std::unordered_map<read_aln_t, bool, read_aln_t::hash_t>> reads;
  };

  // Each node's coverage info
  std::unordered_map<node_t*, aln_info_t> m_coverage_info;

  // History of prior nodes
  typedef std::deque<edge_t*> history_t;

  // A point is an in progress trace.  It is attached to some node, and it
  // holds a seqset range state
  struct point_t {
    point_t(const seqset_range& _range, node_t* _node, const history_t& _history)
        : range(_range), node(_node), history(_history) {}
    bool operator<(const point_t& rhs) const;
    std::string as_string() const;

    seqset_range range;  // The current range
    node_t* node;        // The current node
    history_t history;   // Paths of the range
  };

  // Nodes reachable from a point, and minimum/maximum distance.
  typedef std::map<node_t*, int32_t> reachable_t;

  // Split a reference node at position
  void split_ref(uint32_t pos);

  // Add in the coverage data for a read
  void add_read(const readmap& rm, const point_t& p, const reachable_t& reachable, uint32_t off,
                uint32_t read_id);

 public:
  // Construct a vargraph for a contig
  //vargraph(const dna_slice& contig, size_t min_pair = 100, size_t max_pair = 1000);
  vargraph(const dna_sequence& contig, size_t min_pair = 100, size_t max_pair = 1000);

  // Clean up
  ~vargraph();

  // Load in a variant
  // Here we allow empty sequence for deletions
  void add_variant(uint32_t start, uint32_t end, const dna_sequence& seq);

  // Trace everything
  void trace(const seqset& ss, const readmap& rm, uint32_t start, uint32_t end);

  // Get all the entries
  const std::multimap<uint32_t, node_t*>& get_nodes() { return m_gnodes; }

  // Get the edges 
  const std::set<edge_t*>& get_edges() { return m_edges; }
};
