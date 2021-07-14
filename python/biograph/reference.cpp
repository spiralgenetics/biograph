#include <pybind11/pybind11.h>

#include "modules/bio_base/seq_position.h"
#include "python/biograph/reference.h"

////////
// This is part of the public python API.
//
// Note that all ranges are 0-based internally
////////

// Make a new reference_range. If use_exact_loci is false, attempt to move the
// start or end outside of any containing N region. Throws if start and end
// are not within the same supercontig.
reference_range reference_wrapper::make_range(const std::string& supercontig_name,
                                              unsigned long start, unsigned long end,
                                              bool use_exact_loci) const {
  std::pair<size_t, size_t> flat_start_end =
      m_reference->flatten_range(supercontig_name, start, end, use_exact_loci);

  return reference_range(m_reference, flat_start_end.first, flat_start_end.second);
}

// The number of bases in the whole reference_assembly
size_t reference_wrapper::size() const { return m_reference->size(); }

// Return a py_list of supercontigs in order.
py_list reference_wrapper::supercontigs() const {
  std::vector<std::string> supercontigs = m_reference->get_assembly().get_supercontig_order();

  py_list list;
  for (const auto& supercontig : supercontigs) {
    // This may need to be moved.. why do we go all through every time?
    // Seems wasteful
    size_t colon = supercontig.find_last_of(":");

    list.append(
        get_supercontig(supercontig.substr(0, colon), stoul(supercontig.substr(colon + 1))));
  }

  return list;
}

// Return a py_list of scaffolds
py_list reference_wrapper::scaffolds() const {
  std::set<scaffold> scaffolds = m_reference->get_assembly().scaffolds;

  py_list list;
  for (const auto& s : scaffolds) {
    list.append(s.name);
  }
  return list;
}

// Return the containing supercontig at a 0-based position as a reference_range
reference_range reference_wrapper::get_supercontig(std::string supercontig_name, size_t pos) const {
  const reference_assembly k = m_reference->get_assembly();

  supercontig contig = k.get_supercontig(
      m_reference->flatten(seq_position(k.get_scaffold(supercontig_name).index, pos)));
  // Hack.2
  // return make_range(contig.scaffold_name, contig.offset + 1, contig.offset +
  // contig.len, false);
  return make_range(contig.scaffold_name, contig.offset, contig.offset + contig.len - 1, false);
}

py_list reference_wrapper::find_ranges(const std::string& supercontig_name, size_t start,
                                       size_t end) const {
  py_list results;
  std::vector<std::string> supercontigs = m_reference->get_assembly().get_supercontig_order();
  for (const auto& supercontig : supercontigs) {
    size_t colon = supercontig.find_last_of(":");
    if (supercontig.substr(0, colon) != supercontig_name) {
      continue;
    }
    auto super_contig_range =
        get_supercontig(supercontig.substr(0, colon), stoul(supercontig.substr(colon + 1)));

    if (super_contig_range.get_start() > end || super_contig_range.get_end() < start) {
      continue;
    }
    results.append(make_range(supercontig_name, std::max(super_contig_range.get_start(), start),
                              std::min(super_contig_range.get_end(), end), false));
  }

  return results;
}

py_dict reference_wrapper::scaffold_lens() const {
  std::set<scaffold> scaffolds = m_reference->get_assembly().scaffolds;

  py_dict ret;
  for (const auto& s : scaffolds) {
    ret[s.name.c_str()] = s.len;
  }
  return ret;
}

// String representation in [inclusive-exclusive) notation
std::string reference_wrapper::str() const { return path(); }

// Python object repr()
std::string reference_wrapper::repr() const {
  return std::string("<biograph.Reference ") + str() + std::string(">");
}

bwt_wrapper reference_wrapper::find_sequence(const dna_sequence& seq) const {
  // if(seq == dna_sequence("")) {
  //     throw io_exception("Cannot call find() on an empty sequence");
  // }
  return bwt_wrapper(m_reference, seq);
}

bwt_wrapper reference_wrapper::find(const std::string& seq_str) const {
  return find_sequence(dna_sequence(seq_str));
}

// A reference_range represents a region of the reference in 0-based
// coordinates. A given range cannot cross supercontig boundaries.
reference_range::reference_range(std::shared_ptr<reference> the_reference, size_t flat_start,
                                 size_t flat_end)
    : m_reference(the_reference), m_flat_ref_start(flat_start), m_flat_ref_end(flat_end) {}

// String representation in [inclusive-exclusive) notation
std::string reference_range::str() const {
  return printstring("%s:[%lu-%lu)", get_scaffold().c_str(), get_start(), get_end());
}

// Python object repr()
std::string reference_range::repr() const {
  return std::string("<biograph.ReferenceRange ") + str() + std::string(">");
}

// The name of the scaffold containing this range.
std::string reference_range::get_scaffold() const {
  return m_reference->get_assembly()
      .scaffold_order[m_reference->get_seq_position(m_flat_ref_start).scaffold_id];
}

// First absolute (*seqpos*) position in this range.
unsigned long reference_range::get_start() const {
  return m_reference->get_seq_position(m_flat_ref_start).position;
}

// Last absolute (*seqpos*) position in this range.
unsigned long reference_range::get_end() const {
  return m_reference->get_seq_position(m_flat_ref_end).position;
}

// The sequence represented by this range.
dna_sequence reference_range::sequence() const {
  return dna_sequence(m_reference->get_dna(m_flat_ref_start), m_reference->get_dna(m_flat_ref_end));
}

// Number of bases in this range.
size_t reference_range::size() const { return m_flat_ref_end - m_flat_ref_start; }

// Return a reference_range for the given match
// throw (raise a python RuntimeError) if the bwt is not valid
reference_range bwt_wrapper::get_match(size_t which) const {
  if (!valid()) {
    throw io_exception("Called get_match() on an invalid bwt object");
  }
  if (which >= matches()) {
    throw io_exception(printstring("Called get_match() on invalid match %lu", which));
  }
  size_t start = m_bwt.get_match(which);
  return reference_range(m_reference, start, start + m_query.size());
}

// String representation in [inclusive-exclusive) notation
std::string bwt_wrapper::str() const { return ref_path() + "{" + query().as_string() + "}"; }

// Python object repr()
std::string bwt_wrapper::repr() const {
  return std::string("<biograph.ReferenceContext ") + str() + std::string(">");
}

// pass through all other bwt calls
bool bwt_wrapper::valid() const { return m_bwt.valid(); }

size_t bwt_wrapper::begin() const { return m_bwt.begin(); }

size_t bwt_wrapper::end() const { return m_bwt.end(); }

size_t bwt_wrapper::matches() const { return m_bwt.matches(); }

using namespace pybind11;
void bind_reference(module& m) {
  class_<reference_wrapper>(
      m, "Reference",
      "Open a reference genome directory and return a Reference object.\n"
      "\n"
      "The Reference object can then be quickly searched and used for "
      "analysis.\n"
      "\n"
      "Raises a RuntimeError if the reference cannot be opened.\n"
      "\n"
      "Example:\n"
      "    >>> from biograph import Reference\n"
      "    >>> my_ref = Reference('datasets/lambdaToyData/benchmark/ref_lambda')\n"
      "    >>> human_ref = Reference('/reference/human_g1k_v37')\n"
      "")
      .def(init<const std::string&>())
      .def("__str__", &reference_wrapper::str, "A string representation of the Reference object.")
      .def("__repr__", &reference_wrapper::repr,
           "The python representation of the Reference object.")
      .def("make_range", &reference_wrapper::make_range, arg("scaffold_name"), arg("start"),
           arg("end"), arg("use_exact_loci") = true, R"DOC(
Return a ReferenceRange object for the given region of the reference.

The scaffold name is a string. Sequences start with position 0.

Start is inclusive, end is exclusive.

If **use_exact_loci** is True and either end of the region specifies
an N base, a *RuntimeError* is raised. If **use_exact_loci** is False
and one end of the specified region specifies an N base, it will
be automatically moved to the first base outside the N region.

A *RuntimeError* is always raised if both ends of the region specify N bases.

A RuntimeError is also raised if the region is outside the given
scaffold, or if the scaffold is not present in the reference.

Example:
    >>> ref_range = human_ref.make_range('2', 0, 200000, True)
    Traceback (most recent call last):
      ...
    RuntimeError: 2:0 is in a part of the reference (offset 27731782) with no sequence data (only N bases)

    >>> ref_range = human_ref.make_range('2', 0, 200000, False)
    >>> print(ref_range.start)
    10000
)DOC")
      .def_property_readonly("size", &reference_wrapper::size,
                             "The total number of bases in the reference.")
      .def_property_readonly("scaffolds", &reference_wrapper::scaffolds,
                             "A list of all sequences in the reference. Sequences "
                             "include regions\n"
                             "of contiguous DNA, typically separated by regions of "
                             "unknown (N) bases.\n")
      .def_property_readonly("chromosomes", &reference_wrapper::scaffolds,
                             "A synonym for scaffolds.\n")
      .def_property_readonly("supercontigs", &reference_wrapper::supercontigs,
                             "A list of reference_range objects that include all supercontigs\n"
                             "(contiguous DNA regions) in the reference.")
      .def_property_readonly("scaffold_lens", &reference_wrapper::scaffold_lens,
                             "A dict of chromosomes and their lengths\n")
      .def("get_supercontig", &reference_wrapper::get_supercontig,
           arg("scaffold_name"), arg("position"),
           "Return a reference_range object representing the supercontig that "
           "contains\n"
           "the given position.\n")
      .def("find_ranges", &reference_wrapper::find_ranges, arg("scaffold_name"),
           arg("start"), arg("end"),
           "Find all contigs in the given region. Returns a list of "
           "reference_ranges\n")
      .def("find", &reference_wrapper::find, arg("seq"))
      .def("find", &reference_wrapper::find_sequence, arg("seq"),
           "find(sequence)\n"
           "\n"
           "Search the reference for the given sequence or string. Returns a\n"
           "ReferenceContext object.\n"
           "\n"
           "If the sequence is not found in the reference, "
           "reference_context.valid\n"
           "will be False.\n"
           "\n"
           "Example:\n"
           "   >>> my_ref.find('A')\n"
           "   <biograph.ReferenceContext datasets/lambdaToyData/benchmark/ref_lambda/{A}>\n"
           "   >>> my_ref.find(biograph.Sequence('A'))           # equivalent\n"
           "   <biograph.ReferenceContext datasets/lambdaToyData/benchmark/ref_lambda/{A}>\n"
           "   >>> my_ref.find('ACATATTACGGGGT').valid\n"
           "   False\n"
           "   >>> my_ref.find('XYZ')\n"
           "   Traceback (most recent call last):\n"
           "      ...\n"
           "   RuntimeError: Failed conversion of dna_base, c = 'X'\n");

  class_<reference_range>(
      m, "ReferenceRange",
      "An object representing a range of bases within a reference. NOTE: reference locations\n"
      "are consistent with Python's zero-based half-open indexing.\n"
      "")
      .def("__str__", &reference_range::str,
           "A string representation of the specified range. [ indicates that \n"
           "the start is inclusive, while ) indicates that the end is "
           "exclusive.\n"
           "\n"
           "Example:\n"
           "    >>> ref_range = human_ref.make_range('22', 44682664, 44682668, True)\n"
           "    >>> print(ref_range)\n"
           "    22:[44682664-44682668)")
      .def("__repr__", &reference_range::repr,
           "The python representation of this ReferenceRange object.\n"
           "\n"
           "Example:\n"
           "\n"
           "    >>> ref_range = human_ref.make_range('22', 44682664, 44682668, True)\n"
           "    >>> print(repr(ref_range))\n"
           "    <biograph.ReferenceRange 22:[44682664-44682668)>\n")
      .def_property_readonly("scaffold", &reference_range::get_scaffold,
                             "The name of the scaffold used for this range.\n")
      .def_property_readonly("chromosome", &reference_range::get_scaffold,
                             "A synonym for scaffold.\n")
      .def_property_readonly("start", &reference_range::get_start,
                             "The starting position of this range (zero-based, inclusive)\n")
      .def_property_readonly("end", &reference_range::get_end,
                             "The end position of this range (zero-based, exclusive)\n")
      .def_property_readonly(
          "sequence", &reference_range::sequence,
          "A Sequence object representing the nucleotide sequence\n"
          "for this range.\n"
          "\n"
          "Example:\n"
          "    >>> ref_range = human_ref.make_range('22', 44682664, 44682668, True)\n"
          "    >>> print(ref_range.sequence)\n"
          "    GGGC\n")
      .def_property_readonly("size", &reference_range::size,
                             "The number of bases in this range.\n");

  class_<bwt_wrapper>(m, "ReferenceContext",
                      "A ReferenceContext object representing a set of "
                      "reference locations that\n"
                      "match a query, returned by Reference.find().\n"
                      "")
      .def("__str__", &bwt_wrapper::str,
           "A string representation of the ReferenceContext object.\n")
      .def("__repr__", &bwt_wrapper::repr,
           "The python representation of the ReferenceContext object.\n")
      .def_property_readonly("matches", &bwt_wrapper::matches,
                             "The number of times that query matches the reference.\n")
      .def_property_readonly("valid", &bwt_wrapper::valid,
                             "True if query matches the reference in at least one location,\n"
                             "otherwise False.\n")
      .def_property_readonly("start", &bwt_wrapper::begin, "The first matching reference entry.\n")
      .def_property_readonly("end", &bwt_wrapper::end, "The last matching reference entry.\n")
      .def("get_match", &bwt_wrapper::get_match,
           "Return a ReferenceContext for the given matching entry.\n")
      .def_property_readonly("query", &bwt_wrapper::query,
                             "The query sequence used to generate this ReferenceContext.\n");
}
