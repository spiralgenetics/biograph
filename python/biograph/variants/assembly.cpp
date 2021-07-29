#include "python/biograph/variants/assembly.h"

#include "modules/bio_base/readmap.h"
#include "modules/variants/assemble.h"
#include "python/common.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

using namespace pybind11;
using namespace variants;

bool read_coverage_contains(const read_coverage_t& c, const read_coverage_read_t& r) {
  const auto& reads = c.reads();
  auto eq = std::equal_range(reads.begin(), reads.end(), r, read_coverage_read_order());
  if (eq.first == eq.second) {
    return false;
  }

  std::unordered_set<uint32_t> ids_needed(r.read_ids.begin(), r.read_ids.end());
  for (auto it = eq.first; it != eq.second; ++it) {
    for (uint32_t read_id : it->read_ids) {
      ids_needed.erase(read_id);
      if (ids_needed.empty()) {
        return true;
      }
    }
  }

  return false;
}

std::string read_coverage_read_repr(read_coverage_read_t& rd) {
  std::string read_ids;
  if (rd.read_ids.size() == 1) {
    for (uint32_t read_id : rd.read_ids) {
      read_ids = std::to_string(read_id);
    }
  } else {
    for (uint32_t read_id : rd.read_ids) {
      if (read_ids.empty()) {
        read_ids += "[";
      } else {
        read_ids += ",";
      }
      read_ids += std::to_string(read_id);
    }
    if (read_ids.empty()) {
      read_ids = "(none)";
    } else {
      read_ids += "]";
    }
  }
  return "ReadCoverageRead(" + std::to_string(rd.offset) + ", " + read_ids + ", " +
         std::to_string(rd.read_len) + ")";
}

std::vector<read_coverage_read_t>::const_iterator read_coverage_begin(
    const read_coverage_t& reads) {
  return reads.reads().begin();
}

std::vector<read_coverage_read_t>::const_iterator read_coverage_end(const read_coverage_t& reads) {
  return reads.reads().end();
}

std::vector<int> read_coverage_calc_depths(const read_coverage_t& reads, bool include_fwd,
                                           bool include_rev, bool interbase,
                                           const std::shared_ptr<readmap>& rm) {
  return reads.calc_depths(include_fwd, include_rev, interbase, rm.get());
}

list read_id_set_expand_to_list(const read_id_set& reads) {
  list result;
  for (uint32_t read_id : reads) {
    result.append(read_id);
  }
  return result;
}

void bind_assembly(module& m) {
  class_<read_id_set>(m, "ReadIdSet", "Contains a set of read ids")
      .def(init<>())
      .def(init<>([](object input_reads) {
        read_id_set reads;
        for (auto elem : input_reads) {
          reads.insert(cast<int>(elem));
        }
        return reads;
      }))
      .def("add", (void (read_id_set::*)(uint32_t)) & read_id_set::insert)
      .def("expand_to_list", read_id_set_expand_to_list)
      .def("__len__", &read_id_set::size)
      .def("__str__", str_from_ostream<read_id_set>)
      .def(
          "__iter__",
          [](const read_id_set& ids) {
            return make_iterator<return_value_policy::copy>(ids.begin(), ids.end());
          },
          keep_alive<0, 1>())
      .def(self + self)  // TODO(nils): Deprecate this in favor of bitwise-or operator
      .def(self == self)
      .def(self | self)
      .def(self - self)
      .def(self & self)
      .def(self |= self)
      .def(self &= self)
      .def(self -= self);

  class_<big_read_id_set>(
      m, "BigReadIdSet",
      "Contains a set of read ids, optimized for use when many read ids are present.")
      .def(init<>())
      .def("add", (void (big_read_id_set::*)(uint32_t)) & big_read_id_set::insert)
      .def("__len__", &big_read_id_set::size)
      .def("__str__", str_from_ostream<big_read_id_set>)
      .def(
          "__iter__",
          [](const big_read_id_set& ids) {
            return make_iterator<return_value_policy::copy>(ids.begin(), ids.end());
          },
          keep_alive<0, 1>())
      .def("to_read_id_set", [](const big_read_id_set& ids) { return read_id_set(ids); })
      .def(self | read_id_set())
      .def(self - read_id_set())
      .def(self & read_id_set())
      .def(self |= read_id_set())
      .def(self &= read_id_set())
      .def(self -= read_id_set());

  class_<edge_coverage_t>(m, "EdgeCoverage",
                          "Contains information on edge coverage for an assembly.")
      .def(init<>())
      .def("__str__", str_from_ostream<edge_coverage_t>)
      .def("__repr__", str_from_ostream<edge_coverage_t>)
      .def_readonly("variant_start", &edge_coverage_t::variant_start,
                    "Coverage for the assembly's left anchor")
      .def_readonly("variant_end", &edge_coverage_t::variant_end,
                    "Coverage for the assembly's right anchor")
      .def_readonly("reference_start", &edge_coverage_t::reference_start,
                    "Coverage for the reference at the assembly's left anchor")
      .def_readonly("reference_end", &edge_coverage_t::reference_end,
                    "Coverage for the reference at the assembly's right anchor");

  class_<align_count_t>(m, "AlignCount",
                          "Contains information on alignment counts for overlapping assemblies.")
      .def(init<>())
      .def("__str__", str_from_ostream<align_count_t>)
      .def_readonly("local_aligned_bases", &align_count_t::local_aligned_bases,
                    "Sum of lengths of first aligments of reads in this assembly")
      .def_readonly("local_read_lens", &align_count_t::local_read_lens,
                    "Sum of read lengths of all reads with alignments in this assembly")
      .def_readonly("tot_aligned_bases", &align_count_t::tot_aligned_bases,
                    "Sum of all alignments for this read overlapping this assembly");

  class_<read_coverage_read_t>(m, "ReadCoverageRead",
                               "A single read aligned to a specific position in an assembly")
      .def(init<aoffset_t, uint32_t, int>())
      .def_readwrite("offset", &read_coverage_read_t::offset)
      .def_readwrite("read_ids", &read_coverage_read_t::read_ids)
      .def_readwrite("read_len", &read_coverage_read_t::read_len)
      .def("__repr__", read_coverage_read_repr)
      .def("__str__", read_coverage_read_repr)
      .def(self == self);

  class_<read_coverage_t>(m, "ReadCoverage",
                          "Contains information on what reads provide coverage for an assembly.  "
                          "This is a set of (read id, offset of beginning of read from beginning "
                          "of assembly).")
      .def(init<>())
      .def(init([](aoffset_t assembly_len, object input) {
        read_coverage_set new_reads;
        for (auto item : input) {
          new_reads.insert(cast<const read_coverage_read_t&>(item));
        }
        return new_reads.build_and_clear(assembly_len);
      }))
      .def(
          "__iter__",
          [](const read_coverage_t& cov) {
            const auto& reads = cov.reads();
            return make_iterator<return_value_policy::copy>(reads.begin(), reads.end());
          },
          keep_alive<0, 1>())
      .def("__contains__", read_coverage_contains)
      .def("__len__", [](const read_coverage_t& cov) { return cov.reads().size(); })
      .def(self == self)
      .def(self != self)
      .def("assembly_len", &read_coverage_t::assembly_len)
      .def("calc_depths", read_coverage_calc_depths, arg("include_fwd") = true,
           arg("include_rev") = true, arg("interbase") = true,
           arg("readmap") = std::shared_ptr<readmap>())
      .def("get_reads_spanning_offset", &read_coverage_t::get_reads_spanning_offset)
      .def("get_overlaps", &read_coverage_t::get_overlaps)
      .def("get_overlap_min_max", &read_coverage_t::get_overlap_min_max)
      .def("get_and_adjust_reads_spanning_offset",
           &read_coverage_t::get_and_adjust_reads_spanning_offset)
      .def("union_with", &read_coverage_t::union_with)
      .def("intersection_with", &read_coverage_t::intersection_with)
      .def("get_max_flank", &read_coverage_t::get_max_flank)
      .def("get_tot_read_count", &read_coverage_t::get_tot_read_count)
      .def("all_read_ids", &read_coverage_t::all_read_ids)
      .def(self | self)
      .def(self & self)
      .def(self - self)
      .def(self |= self)
      .def(self &= self)
      .def(self -= self)
      .def(self & read_id_set())
      .def(self - read_id_set())
      .def(self &= read_id_set())
      .def(self -= read_id_set());

  class_<assembly, assembly_ptr>(m, "Assembly", "Contains a sequence anchored to reference\n",
                                 dynamic_attr())
      .def(init([](optional_aoffset left_offset, optional_aoffset right_offset,
                   const dna_sequence& seq, size_t assembly_id) -> assembly_ptr {
        return assembly_ptr(make_unique<assembly>(left_offset, right_offset, seq, assembly_id));
      }))
      .def_readwrite("assembly_id", &assembly::assembly_id)
      .def_readwrite("left_offset", &assembly::left_offset,
                     "0-based offset to reference the left side of this "
                     "assembly is anchored to")
      .def_readwrite("right_offset", &assembly::right_offset,
                     "0-based offset to reference the right side of this "
                     "assembly is anchored to")
      .def_readwrite("seq", &assembly::seq, "Sequence of bases")
      .def_readwrite("edge_coverage", &assembly::edge_coverage)
      .def_readwrite("read_coverage", &assembly::read_coverage)
      .def_readwrite("pair_read_coverage", &assembly::pair_read_coverage)
      .def_readwrite("align_count", &assembly::align_count)
      .def_readwrite("matches_reference", &assembly::matches_reference,
                     "True if this assembly matches reference entirely")
      .def_readwrite("min_overlap", &assembly::min_overlap,
                     "Minimum overlap seen along this assembly path")
      .def_readwrite("phase_ids", &assembly::phase_ids)
      .def_property(
          "tags",
          [](const assembly& a) {
            std::vector<std::string> tags(a.tags.begin(), a.tags.end());
            return tags;
          },
          [](assembly& a, std::vector<std::string> new_tags) {
            a.tags.clear();
            a.tags.insert(new_tags.begin(), new_tags.end());
          })
      .def_property(
          "generated_by", [](const assembly& a) { return a.tags.to_string_short(); },
          [](assembly& a, std::string new_gen_by) {
            a.tags.clear();
            a.tags.insert(new_gen_by);
          },
          "Source of this assembly; usually 'PUSH' for push tracer or 'POP' for pop tracer")
      .def_readwrite(
          "bypass_coverage", &assembly::bypass_coverage,
          "If true, do not generate coverage for this assembly or trace paths through it.")
      .def_readwrite("read_cov_max_paths", &assembly::read_cov_max_paths)
      .def(
          "reverse_in_place",
          [](assembly* a, aoffset_t ref_end_pos, std::shared_ptr<readmap> rm) {
            reverse_assembly_in_place(a, rm.get(), ref_end_pos);
          },
          arg("ref_end_pos"), arg("readmap") = nullptr,
          "Returns a backwards version of this assembly, where an offset of 0 is mapped to the "
          "given offset, and the given offset is mapped to an offset of 0.")
      .def("__str__", str_from_ostream<assembly>)
      .def("__lt__",
           [](const assembly& a, const assembly& b) { return canon_assembly_order()(a, b); })
      .def("__ge__",
           [](const assembly& a, const assembly& b) { return !canon_assembly_order()(a, b); });
}
