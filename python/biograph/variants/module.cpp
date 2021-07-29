#include "python/biograph/variants/module.h"

#include <pybind11/pybind11.h>

#include "python/biograph/variants/add_ref.h"
#include "python/biograph/variants/align_reads.h"
#include "python/biograph/variants/apply_edges.h"
#include "python/biograph/variants/apply_graph.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/dedup_cov_reads.h"
#include "python/biograph/variants/discover.h"
#include "python/biograph/variants/filter_dup_align.h"
#include "python/biograph/variants/graph_discover.h"
#include "python/biograph/variants/limit_alleles.h"
#include "python/biograph/variants/pair_cov.h"
#include "python/biograph/variants/pair_edge_cov.h"
#include "python/biograph/variants/phase.h"
#include "python/biograph/variants/place_pair_cov.h"
#include "python/biograph/variants/read_cov.h"
#include "python/biograph/variants/trim_ref.h"
#include "python/common.h"

using namespace variants;

using namespace pybind11;

void bind_variants_module(module& m) {
  bind_assembly(m);
  bind_pair_edge_cov(m);
  bind_read_cov(m);
  bind_dedup_cov_reads(m);
  bind_trim_ref(m);
  bind_pair_cov(m);
  bind_add_ref(m);
  bind_discover(m);
  bind_apply_edges(m);
  bind_phases(m);
  bind_limit_alleles(m);
  bind_filter_dup_align(m);
  bind_align_reads(m);
  bind_graph_discover(m);
  bind_place_pair_cov(m);
  bind_apply_graph(m);
}
