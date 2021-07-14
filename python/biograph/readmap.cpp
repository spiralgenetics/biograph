#include "python/biograph/readmap.h"
#include <pybind11/operators.h>
#include "modules/bio_base/dna_sequence.h"
#include "python/biograph/seqset.h"
#include "python/common.h"

using namespace pybind11;

namespace {

std::string readmap_repr(const readmap& rm) {
  return printstring("<Readmap %p: %s>", &rm, rm.path().c_str());
}

std::string readmap_read_repr(const readmap::read& read) {
  return printstring("<ReadmapRead id=%d>", read.get_read_id());
}

size_t readmap_read_hash(const readmap::read& r) { return std::hash<uint32_t>()(r.get_read_id()); }

std::vector<int> readmap_get_approx_seq_coverage(const readmap& rm, dna_sequence seq) {
  return rm.approx_coverage(dna_slice(seq.begin(), seq.end()));
}

}  // namespace

void bind_readmap(module& m) {
  class_<readmap, std::shared_ptr<readmap>>(m, "Readmap",
                  R"DOC(
Readmap
Load a Readmap. Returns a new Readmap object loaded in memory and
ready to query.

Calling Readmap with an empty string returns a dummy Readmap
representing all possible reads.
This rarely needs to be called directly since the BioGraph object is
a combination of the Seqset and Readmap and holds the main operations.

Raises a RuntimeError if the file cannot be opened.
)DOC")
      .def("__repr__", readmap_repr)
      .def("__str__", readmap_repr)
      .def_property_readonly("seqset", &readmap::get_seqset,
                             "The seqset this Readmap is assoiated with.")
      .def("max_read_len", &readmap::max_read_len, "Maximum read length present in this readmap")
      .def("min_read_len", &readmap::min_read_len, "Minimum read length present in this readmap")
      .def("size", &readmap::size,
           "How many read IDs are present in this readmap.  (Reads are tracked in both forward "
           "and "
           "reverse complement form, so this will be twice the number of actual reads)")
      .def("get_read_by_id", &readmap::get_read_by_id, arg("id"), "Lookup a read by read id")
      .def("get_prefix_reads", &readmap::get_prefix_reads, arg("entry"), arg("min_read_len") = 0,
           "Iterate through reads that are a prefix of a given SeqsetEntry")
      .def("get_reads_containing", &readmap::get_reads_containing, arg("entry"),
           "Iterate through reads that contain a given SeqsetEntry somewhere inside the read. "
           "Returns a tuple of (offset, ReadmapRead) where offset is the number of bases before "
           "the given SeqsetEntry in the read.")
      .def("get_approx_seq_coverage", readmap_get_approx_seq_coverage, arg("seq"),
           "Returns the coverage of a sequence. Sequence must be longer than read length.  This "
           "may return the wrong value in some cases where the read length is shorter than the "
           "seqset entry length.")
      .def("get_read_count", &readmap::get_read_count, "Return the number of reads in this readmap")
      .def("get_num_bases", &readmap::get_num_bases, "Return the number of bases in this readmap")
      .def("get_pair_stats", &readmap::get_pair_stats,
           "Return the number of paired/unpaired reads/bases");

  class_<readmap::pair_stats>(m, "ReadmapPairStats",
                              R"DOC(
ReadmapPairStats
Statics on pairing data within a readmap
)DOC")
      .def_readwrite("paired_reads", &readmap::pair_stats::paired_reads,
                     "Number of paired reads present")
      .def_readwrite("paired_bases", &readmap::pair_stats::paired_bases,
                     "Number of bases present in paired reads")
      .def_readwrite("unpaired_reads", &readmap::pair_stats::unpaired_reads,
                     "Number of unpaired reads present")
      .def_readwrite("unpaired_bases", &readmap::pair_stats::unpaired_bases,
                     "Number of bases present in unpaired reads");

  class_<readmap::read_iterator_range>(m, "_ReadmapReadRange",
                                       R"DOC(
A _ReadmapReadRange iterator representing a range of ReadmapRead entries.
)DOC")
      .def("__iter__", [](const readmap::read_iterator_range& reads) -> object {
          return pybind11::make_iterator<return_value_policy::copy>(reads.begin(), reads.end());
      });

  class_<readmap::containing_read_iterator_range>(m, "_ReadmapContainingReadRange",
                                                  R"DOC(
A _ReadmapContainingReadRange iterator representing a range of ReadmapRead entries.
)DOC")
      .def("__iter__", [](const readmap::containing_read_iterator_range& reads) -> object {
        return pybind11::make_iterator<return_value_policy::copy>(reads.begin(), reads.end());
      });

  class_<readmap::read>(m, "ReadmapRead",
                        R"DOC(
A ReadmapRead object holds properties of a single read.
)DOC")
      .def(self == self)
      .def(self != self)
      .def("__len__", &readmap::read::size, "Length of this read")
      .def("__repr__", readmap_read_repr)
      .def("__str__", readmap_read_repr)
      .def("__hash__", readmap_read_hash)
      .def("has_mate", &readmap::read::has_mate, "Returns true if this read is part of a pair")
      .def("get_mate", &readmap::read::get_mate, "Returns the mate associated with this read")
      .def("get_rev_comp", &readmap::read::get_rev_comp,
           "Returns the representation of this read that's the reverse complement of this read")
      .def("is_original_orientation", &readmap::read::is_original_orientation,
           "Is this read the direction we originally saw it?")
      .def("get_seqset_entry", &readmap::read::get_seqset_entry,
           "Seqset entry associated with this read")
      .def("get_read_id", &readmap::read::get_read_id, "Returns the numeric read id of this read");
}
