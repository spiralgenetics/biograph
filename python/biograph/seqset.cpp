#include "python/biograph/seqset.h"
#include "python/biograph/readmap.h"
#include "python/biograph/reference.h"
#include "python/common.h"

#include <pybind11/functional.h>
#include <pybind11/operators.h>

using namespace pybind11;

namespace {

size_t do_seqset_range_hash(const seqset_range& r) { return seqset_range_hash()(r); }

std::string seqset_range_repr(const seqset_range& r) {
  return "<SeqsetEntry " + std::to_string(r.begin()) + "-" + std::to_string(r.end()) + ": " +
         r.sequence().as_string() + ">";
}

template <typename... Arg>
std::function<object(Arg...)> valid_or_none_f(std::function<seqset_range(Arg... arg)> f) {
  return [f](Arg... arg) -> object {
    seqset_range r = f(arg...);
    if (r.valid()) {
      return cast(r);
    } else {
      return none();
    }
  };
}

template <typename T, typename... Arg>
std::function<object(const T&, Arg...)> valid_or_none(seqset_range (T::*f)(Arg... args) const) {
  return [f](const T& obj, Arg&&... arg) -> object {
    seqset_range r = (obj.*f)(std::forward<Arg>(arg)...);
    if (r.valid()) {
      return cast(r);
    } else {
      return none();
    }
  };
}

}  // namespace

void bind_seqset(module& m) {
  class_<seqset, std::shared_ptr<seqset>>(m, "Seqset",
                                          R"DOC(
Seqset
Load a Seqset. Returns a new Seqset object loaded in memory and ready
to query.

Seqsets contain the read overlap graph, which can be used for 
arbitrary kmer queries. For depth of coverage queries, you should 
load the Readmap that was generated for the Seqset.

This rarely needs to be called directly since the BioGraph object is 
a combination of the Seqset and Readmap and holds the main operations.

Raises a RuntimeError if the file cannot be opened.

Example:

    from biograph import Seqset
    my_sample = Seqset('/datasets/na12878.seqset')

See also: Readmap
)DOC")
      .def("size", &seqset::size,
           R"DOC(
Total number of entries in the Seqset.

Example:
    print my_sample.size  # 2106729563
)DOC")
      .def("find",
           valid_or_none((seqset_range(seqset::*)(const dna_sequence&) const)(&seqset::find)),
           arg("self"),  // arg("seq"),
           R"DOC(
Search the Seqset for this sequence and return a new SeqsetEntry 
object representing the search results. The query can be a sequence 
object or string of nucleotides. If the query is not found in the 
Seqset, the resulting SeqsetEntry object's 'valid' attribute will 
be False. Otherwise, 'valid' will be True and the matching locations 
can be accessed.

Example:

    entry = my_sample.find('ACGT')
    entry = my_sample.find(sequence('ACGT')) # equivalent
    print entry.valid  # True

    bad_entry = my_sample.find('TTT').next_kmer()
    print bad_entry.valid  # False

See also: SeqsetEntry
)DOC")
      .def("empty_entry", valid_or_none(&seqset::ctx_begin),
           R"DOC(
Returns a SeqsetEntry representing an empty sequence; this is a 
prefix of all entries in the seqset.  This is equivalent to passing 
an empty sequence to find().

)DOC")
      .def("get_entry_by_id", valid_or_none(&seqset::ctx_entry), arg("id"),
           R"DOC(
Returns a new SeqsetEntry representing the single specified 
entry.

Example:

    # Print all maximal entries starting with this sequence
    entry = my_sample.find('CATTTAGGACACCT')
    for i in range(entry.start, entry.end):
        print my_sample.get_entry(i).sequence
)DOC")
      .def("max_sequence_length", &seqset::max_read_len,
           R"DOC(
Returns the maximum length of any sequence present in the seqset.

Example:
  >>> my_seqset = my_biograph.seqset
  >>> max_size = max(len(my_seqset.get_entry_by_id(id)) for id in range(my_seqset.size()))
  >>> max_size == my_seqset.max_sequence_length()
  True
)DOC");

  // Invalid ranges are represented as python's None; all methods
  // returning a SeqsetEntry should use valid_or_none()
  class_<seqset_range>(m, "SeqsetEntry",
                       R"DOC(
A SeqsetEntry object representing a sequence contained in Seqset and 
represented as a [start-end) range.
)DOC")
      .def("sequence", &seqset_range::sequence)
      .def("__len__", &seqset_range::size)
      .def("truncate", valid_or_none(&seqset_range::truncate))
      .def("push_front",
           valid_or_none_f(std::function<seqset_range(const seqset_range&, char)>(
               [](const seqset_range& r, char b) -> seqset_range {
                 return r.push_front(dna_base(b));
               })),
           arg("base"),
           R"DOC(
Return a new entry with the given base pushed to the front of the 
current sequence.
)DOC")
      .def("sequence", &seqset_range::sequence, arg("bases") = std::numeric_limits<int>::max(),
           R"DOC(
A Sequence object representing the nucleotide sequence of this entry.
)DOC")
      .def("get_begin_entry_id", &seqset_range::begin, "The index of this entry's beginning.")
      .def("get_end_entry_id", &seqset_range::end, "The index of this entry's end.")
      .def("push_front", valid_or_none(&seqset_range::push_front), arg("base"))
      .def("push_front_drop",
           valid_or_none_f(std::function<seqset_range(const seqset_range&, char, int)>(
               [](const seqset_range& r, char b, int min_ctx) -> seqset_range {
                 return r.push_front_drop(dna_base(b), min_ctx);
               })),
           arg("base"), arg("min_ctx") = 0,
           R"DOC(
Add a dna base to the front, dropping as many bases as needed to make 
a valid result.
)DOC")
      .def("pop_front", valid_or_none(&seqset_range::pop_front),
           R"DOC(
Return a new entry with a single base removed from the front of the 
current sequence.

Raises a RuntimeError if the sequence is empty.
           )DOC")
      .def("truncate", valid_or_none(&seqset_range::truncate), arg("new_size"),
           R"DOC(,
Return a new entry with bases removed from the end of the current 
sequence.

)DOC")
      .def(self == self)
      .def(self != self)
      .def("__str__", seqset_range_repr)
      .def("__repr__", seqset_range_repr)
      .def("__hash__", do_seqset_range_hash);
}
