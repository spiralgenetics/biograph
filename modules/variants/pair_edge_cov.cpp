#include "modules/variants/pair_edge_cov.h"

namespace variants {

namespace {
constexpr bool k_cov_debug = false;
constexpr bool k_extended_stats = false;
constexpr int k_report_seconds = 300;
}  // namespace

pair_edge_cov::~pair_edge_cov() { flush(); }

void pair_edge_cov::on_advance(aoffset_t offset) {
  if (k_extended_stats) {
    static time_t last_report = 0;
    time_t now = time(0);
    if (last_report == 0) {
      last_report = now - (k_report_seconds / 2);
    } else if (last_report + k_report_seconds < now) {
      last_report = now;
      std::cout << "Pair edge cov at " << offset << ": " << sorted_output_stats(offset) << "\n";
    }
  }
}

void pair_edge_cov::on_assembly_edges(optional_aoffset cur_pos,
                                      const std::vector<assembly_ptr>& left_edges,
                                      const std::vector<assembly_ptr>& inserts,
                                      const std::vector<assembly_ptr>& right_edges) {
  read_coverage_t ref_ending_here;
  for (const auto& a : left_edges) {
    if (!a->matches_reference) {
      continue;
    }

    ref_ending_here = ref_ending_here.union_with(
        a->read_coverage->get_and_adjust_reads_spanning_offset(a->seq.size()));
  }

  read_coverage_t ref_starting_here;
  for (const auto& a : right_edges) {
    if (!a->matches_reference) {
      continue;
    }
    ref_starting_here =
        ref_starting_here.union_with(a->read_coverage->get_and_adjust_reads_spanning_offset(0));
  }

  read_id_set reference_reads = ref_starting_here.intersection_with(ref_ending_here).all_read_ids();

  for (const auto& a : left_edges) {
    add_var_edge_coverage(a.get());
    if (a->matches_reference) {
      continue;
    }
    a->edge_coverage->reference_end = reference_reads;
  }

  for (const auto& a : inserts) {
    add_var_edge_coverage(a.get());
    CHECK(!a->matches_reference) << "Insert says it matches reference: " << *a;
    a->edge_coverage->reference_start = reference_reads;
    a->edge_coverage->reference_end = reference_reads;
  }

  for (const auto& a : right_edges) {
    add_var_edge_coverage(a.get());
    if (a->matches_reference) {
      continue;
    }
    a->edge_coverage->reference_start = reference_reads;
  }
}

void pair_edge_cov::on_assembly(assembly_ptr a) {
  if (!a->pair_read_coverage) {
    throw(io_exception("pair_edge_cov requires pair_cov to be run first"));
  }
  a->edge_coverage.emplace();
  apply_edges_step::on_assembly(std::move(a));
}

pair_edge_cov::pair_edge_cov(const assemble_options& options, pipeline_step_t output)
    : apply_edges_step(std::move(output)), m_options(options) {}

void pair_edge_cov::add_var_edge_coverage(assembly* a) {
  edge_coverage_t& ec = *a->edge_coverage;
  aoffset_t seq_end = a->seq.size();
  for (const auto& cov_entry : a->pair_read_coverage->reads()) {
    aoffset_t start = cov_entry.offset;
    aoffset_t end = cov_entry.offset + cov_entry.read_len;

    if (start < 0) {
      CHECK_GT(end, 0);
      ec.variant_start.insert(cov_entry.read_ids.begin(), cov_entry.read_ids.end());
    }
    if (end > seq_end) {
      CHECK_LT(start, seq_end);
      ec.variant_end.insert(cov_entry.read_ids.begin(), cov_entry.read_ids.end());
    }
    if (start >= 0 && end <= seq_end) {
      ec.interior.insert(cov_entry.read_ids.begin(), cov_entry.read_ids.end());
    }
  }
}

void pair_edge_cov::add_edge_read_ids(const assembly* a, aoffset_t offset, read_id_set* read_ids) {
  for (const auto& cov_entry : a->pair_read_coverage->reads()) {
    aoffset_t start = cov_entry.offset;
    aoffset_t end = cov_entry.offset + cov_entry.read_len;

    if (start >= offset || end <= offset) {
      continue;
    }

    read_ids->insert(cov_entry.read_ids.begin(), cov_entry.read_ids.end());
  }
}

}  // namespace variants
