
#include "vargraph.h"
#include <boost/format.hpp>

void vargraph::edge_t::flatten() {
  paired = 0;
  unpaired = 0;
  for(const auto& kvp : coverage) {
    for(const auto& kvp2 : kvp.second) {
      if (kvp2.second) {
        paired++;
      } else {
        unpaired++;
      }
    }
  }
}

vargraph::read_aln_t vargraph::read_aln_t::make_aln(unsigned& read_size, uint32_t pos) {
  read_aln_t r;
  r.target_end = pos;
  r.query_end = read_size;
  if (read_size > pos) {
    r.target_start = 0;
    r.query_start = read_size - pos;
    read_size -= pos;
  } else {
    r.target_start = pos - read_size;
    r.query_start = 0;
    read_size = 0;
  }
  return r;
}

bool vargraph::read_aln_t::operator==(const read_aln_t& rhs) const {
  return target_start == rhs.target_start && target_end == rhs.target_end &&
         query_start == rhs.query_start && query_end == rhs.query_end;
}

size_t vargraph::read_aln_t::hash_t::operator()(const read_aln_t& x) const {
  return (size_t(x.target_start) << 32) ^ size_t(x.target_end) ^ (size_t(x.query_start) * 99991) ^
         (size_t(x.query_end) * 198091);
}

void vargraph::aln_info_t::flatten(node_t* cur_node) const {
  size_t seq_size = cur_node->seq.size();
  std::vector<int> unpaired_start_count(seq_size + 1);
  std::vector<int> unpaired_end_count(seq_size + 1);
  std::vector<int> paired_start_count(seq_size + 1);
  std::vector<int> paired_end_count(seq_size + 1);
  for (const auto& pos_aln : reads) {
    for (const auto& aln : pos_aln.second) {
      if (aln.second) {// is paired read
        paired_start_count[aln.first.target_start]++;
        paired_end_count[aln.first.target_end]++;
      } else {
        unpaired_start_count[aln.first.target_start]++;
        unpaired_end_count[aln.first.target_end]++;
      }
    }
  }
  int unpaired_cov = 0;
  int paired_cov = 0;
  for (size_t i = 0; i < seq_size; i++) {
    unpaired_cov -= unpaired_end_count[i];
    paired_cov -= paired_end_count[i];
    if (i > 0) {
      cur_node->unpaired.span_cov.push_back(std::min(unpaired_cov, 255));
      cur_node->paired.span_cov.push_back(std::min(paired_cov, 255));
    }
    unpaired_cov += unpaired_start_count[i];
    paired_cov += paired_start_count[i];
    cur_node->unpaired.base_cov.push_back(std::min(unpaired_cov, 255));
    cur_node->paired.base_cov.push_back(std::min(paired_cov, 255));
  }
}

bool vargraph::point_t::operator<(const point_t& rhs) const {
  if (node->start != rhs.node->start) {
    return node->start < rhs.node->start;
  }
  if (node != rhs.node) {
    return node < rhs.node;
  }
  return history < rhs.history;
}

std::string vargraph::point_t::as_string() const {
  std::string str = node->as_string();
  for (const edge_t* e : history) {
    str += "->" + e->upstream->as_string();
  }
  return str;
}

std::string vargraph::node_t::as_string() const {
  std::string ret = boost::str(boost::format("Node [%d, %d) is_ref=%d") % start % end % is_ref);
  if (!is_ref) {
    ret = ret + " seq : " + seq.as_string();
  }
  return ret;
}

std::string vargraph::edge_t::as_string() const {
  return upstream->as_string() + " " + downstream->as_string();
}

void vargraph::split_ref(uint32_t pos) {
  auto it = m_ref_nodes.upper_bound(pos);
  it--;
  node_t* r = it->second;
  if (pos != r->start) {
    node_t* r2 = new node_t;
    r2->is_ref = true;
    r2->start = pos;
    r2->end = r->end;
    r2->seq = r->seq.subseq(pos - r->start, r->end - pos);
    r->end = pos;
    r->seq = r->seq.subseq(0, pos - r->start);
    r2->downstream = std::move(r->downstream);
    r->downstream.clear();
    edge_t* edge = new edge_t(r, r2);
    r->downstream.push_back(edge);
    r2->upstream.push_back(edge);
    m_ref_nodes.emplace(pos, r2);
    m_gnodes.emplace(pos, r2);
    m_edges.emplace(edge);
  }
}

void vargraph::add_read(const readmap& rm, const point_t& pnt, const reachable_t& reachable,
                        uint32_t off, uint32_t read_id) {
  bool has_mate = false;
  // Check for read mate
  if (rm.has_mate(read_id)) {
    // Grab the mate ReadMap Entry
    uint32_t mate_id = rm.get_mate_rc(read_id);
    //SPLOG("Read has a mate %u", mate_id);

    // Now check for mate in reachable nodes
    // node_t : min-distance
    for (auto& reach : reachable) {
      node_t* mate_n = reach.first;
      //SPLOG("Reachable node: %s", mate_n->as_string().c_str());
      uint32_t base_back = reach.second + mate_n->seq.size() + off;
      // Check for mate in the edges
      for (edge_t *e : mate_n->upstream) {
        auto it = e->coverage.find(mate_id);
        if (it != e->coverage.end()) {
          for (auto& cov_read : it->second) {
            uint32_t dist = base_back + cov_read.first;
            //SPLOG("in edge %s we have dist %d + %d = %u", e->as_string().c_str(), base_back,
                  //cov_read.first, dist);
            if (m_min_pair <= dist && dist <= m_max_pair) {
              has_mate = true;
              cov_read.second = true;  // Make it mated!
            }
          }
        }
      } // End edge check

      // Look for the mate node m_coverage_info
      aln_info_t& ci = m_coverage_info[mate_n];
      auto it = ci.reads.find(mate_id);
      if (it != ci.reads.end()) {
        for (auto& cov_read : it->second) {
          // Compute pairing distance!
          uint32_t dist = base_back - cov_read.first.target_start + cov_read.first.query_start;
          //SPLOG("dist of read %d - %d + %d = %u", base_back, cov_read.first.target_start,
                //cov_read.first.query_start, dist);
          //if (dist <= m_max_pair) {
          if (m_min_pair <= dist && dist <= m_max_pair) {
            has_mate = true;
            cov_read.second = true;  // Make it mated!
          }
        }
      } // end of coverage
    } // End of this reachable
  }

  // Step 1: Add it into our coverage maps
  unsigned read_size = pnt.range.size();  // TODO: This is a lie
  // Insert the coverage entry
  m_coverage_info[pnt.node].reads[read_id][read_aln_t::make_aln(read_size, off)] = has_mate;
  // Insert coverage info for the rest of the read as appropriate
  // We don't do this for edges and my main bug case shows that
  // the edge doesn't have the 'pair' problem
  auto it = pnt.history.begin();
  while (it != pnt.history.end() && read_size > 0) {
    edge_t *e = *it;
    e->coverage[read_id][read_size] = has_mate;
    node_t *n = e->upstream;
    //Is it because we're going back to the same node?
    //What if I did a simple - read : read_aln_t
    // And whenver I come across a reachable, I just
    //SPLOG("Setting paired of %u on %s to %s", read_id, n->as_string().c_str(), has_mate ? "true" : "false");
    m_coverage_info[n].reads[read_id][read_aln_t::make_aln(read_size, n->seq.size())] =
        has_mate;
    it++;
  }
}

//vargraph::vargraph(const dna_slice& contig, size_t _min_pair, size_t _max_pair)
vargraph::vargraph(const dna_sequence& contig, size_t _min_pair, size_t _max_pair)
    : m_contig_size(contig.size()), m_contig(contig), m_min_pair(_min_pair), m_max_pair(_max_pair) {
  node_t* n = new node_t;
  n->is_ref = true;
  n->start = 0;
  n->end = contig.size();
  n->backing = contig;
  n->seq = dna_slice(n->backing.begin(), n->backing.end());
  //n->seq = contig; //dna_slice(n->backing.begin(), n->backing.end());
  m_ref_nodes.emplace(0, n);
  m_gnodes.emplace(0, n);
}

vargraph::~vargraph() {
  for (auto& kvp : m_gnodes) {
    delete kvp.second;
  }
  for (edge_t* e : m_edges) {
    delete e;
  }
}

void vargraph::add_variant(uint32_t start, uint32_t end, const dna_sequence& seq) {
  if (start == end && seq.size() == 0) {
    throw std::runtime_error("Adding empty variant!");
  }
  if (start == 0 || end >= m_contig_size) {
    throw std::runtime_error("Variant is not in contig interior.");
  }
  // Split at start and get ref upstream of start
  split_ref(start);
  node_t* upstream_ref = m_ref_nodes.at(start)->upstream[0]->upstream;
  // Split at end and get ref downstream of end
  split_ref(end);
  node_t* downstream_ref = m_ref_nodes.at(end);
  // For deletions, we only add an edge
  if (seq.size() == 0) {
    // Special case for deletions
    edge_t* e = new edge_t(upstream_ref, downstream_ref);
    m_edges.emplace(e);
    upstream_ref->downstream.push_back(e);
    downstream_ref->upstream.push_back(e);
    return;
  }
  // Make a new node containing the ALT sequence
  node_t* n = new node_t;
  m_gnodes.emplace(start, n);
  n->is_ref = false;
  n->start = start;
  n->end = end;
  n->backing = seq;
  n->seq = dna_slice(n->backing.begin(), n->backing.end());
  // Make edges to/from reference
  edge_t *e1 = new edge_t(upstream_ref, n);
  edge_t *e2 = new edge_t(n, downstream_ref);
  // Hook everything up
  m_edges.emplace(e1);
  m_edges.emplace(e2);
  n->upstream.push_back(e1);
  n->downstream.push_back(e2);
  upstream_ref->downstream.push_back(e1);
  downstream_ref->upstream.push_back(e2);
}

void vargraph::trace(const seqset& ss, const readmap& rm, uint32_t start, uint32_t end) {
  // todo tracks our trace points, ordered by 'logical ref position'
  //should reset all coverages to 0, I think..
  // This doesn't actually use 'end' and doesn't trace sub-sequences, either
  // e.g. a single node from 1-100, tracing 40-60 would still get coverage over
  //  all 1-100 bases...
  std::map<point_t, reachable_t> todo;
  // Start on reference at 'start'
  auto it = m_ref_nodes.upper_bound(start);
  it--;
  todo[point_t(ss.ctx_begin(), it->second, {})];
  while (!todo.empty()) {
    // Extract top point
    point_t pnt = todo.begin()->first;
    reachable_t reachable = todo.begin()->second;
    todo.erase(todo.begin());
    //SPLOG("%lu todos", todo.size());

    // Trace till the end
    uint32_t off = 0;
    node_t* cur_node = pnt.node;
    //SPLOG("Processing: %s", cur_node->as_string().c_str());
    // Add myself in as a reachable node to handle within node pairings
    reachable_t self_reachable = reachable;
    self_reachable[cur_node] = -int32_t(cur_node->seq.size()); // I don't understand this.
    //SPLOG("Adding Seq for %s with reachables:", cur_node->as_string().c_str());
    //for (auto it=self_reachable.begin(); it!=self_reachable.end(); ++it) {
      //SPLOG("%s reachable dist %d", it->first->as_string().c_str(), it->second);
    //}
    // Walk over each base
    for (uint32_t i = 0; i < cur_node->seq.size(); i++) {
      // Attach to range
      pnt.range = pnt.range.push_front_drop(dna_base(cur_node->seq[i]).complement());
      off++;
      // If it's maximal, make it into actual reads
      if (pnt.range.is_maximal()) {
        // Get the list of reads
        auto pair = rm.entry_to_index(pnt.range.begin());
        std::stringstream seq;
        seq << pnt.range.sequence();
        // For each read
        for (uint32_t read_id = uint32_t(pair.first); read_id < uint32_t(pair.second); read_id++) {
          //SPLOG("Found a read %u for %s: %s", read_id, pnt.node->as_string().c_str(),
                //seq.str().c_str());
          add_read(rm, pnt, self_reachable, off, read_id);
        }
      }
    }
    //give up on this point.. but this is wrong. So whatever
    if (pnt.range.size() < 5) {
      continue;
    }
    //SPLOG("history of %s", cur_node->as_string().c_str());
    //for (auto it=pnt.history.begin(); it!= pnt.history.end(); ++it) {
      //SPLOG("Current history %s", (*it)->as_string().c_str());
    //}
    auto it = pnt.history.begin();
    //why are we erasing history based on the points' range size?
    ssize_t hist_size = ssize_t(pnt.range.size()) - ssize_t(cur_node->seq.size());
    //SPLOG("hist size of %ld", hist_size);
    while (hist_size > 0 && it != pnt.history.end()) {
      hist_size -= (*it)->upstream->seq.size();
      //SPLOG("hist size now %ld", hist_size);
      it++;
    }
    pnt.history.erase(it, pnt.history.end());
    // Add new point for each new path, auto deduping via set
    for (edge_t* e : cur_node->downstream) {
      // Duplicate existing point
      point_t npnt = pnt;
      // Add this edge to front of history
      npnt.history.push_front(e);
      npnt.node = e->downstream;
      //SPLOG("  Adding / updating point: %s", npnt.as_string().c_str());
      // Insert/find value
      reachable_t& new_reachable = todo[npnt];
      // Add in my current entry
      new_reachable[cur_node] = 0;
      // Add any prior paths (not already there)
      for (const auto& prev : reachable) {
        int32_t new_dist = prev.second + int(cur_node->seq.size());
        //SPLOG("new dist is %d", new_dist);
        if (new_dist < int32_t(m_max_pair)) {
          if (new_reachable.count(prev.first)) {
            new_reachable[prev.first] = std::min(new_reachable[prev.first], new_dist);
          } else {
            new_reachable[prev.first] = new_dist;
          }
        }
      }
    }
  }
  for (const auto& kvp : m_gnodes) {
    m_coverage_info[kvp.second].flatten(kvp.second);
  }
  for (edge_t* e : m_edges) {
    e->flatten();
  }
}

