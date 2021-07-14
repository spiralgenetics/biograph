#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <stdexcept>

#include "modules/bio_base/align_astar.h"
#include "modules/bio_base/dna_sequence.h"
#include "python/biograph/log.h"

using namespace pybind11;

class py_alignment {
 public:
  py_alignment() {}

  size_t size() { return m_vec.size(); }

  object getitem(int i) {
    i += offset;
    if (i < 0 || i >= (int)m_vec.size()) {
      throw std::out_of_range("Invalid index.");
    }
    return make_tuple(m_vec[i].seq_pos, m_vec[i].read_pos);
  }

  object matches() {
    list out;
    int ms_1 = -1;
    int ms_2 = -1;
    int match_len = 0;
    for (size_t i = offset; i < m_vec.size(); i++) {
      if (s1[m_vec[i].seq_pos] == s2[m_vec[i].read_pos]) {
        if (match_len == 0) {
          match_len = 1;
          ms_1 = m_vec[i].seq_pos;
          ms_2 = m_vec[i].read_pos;
        } else {
          match_len++;
        }
      } else {
        if (match_len > 0) {
          out.append(make_tuple(ms_1, ms_2, match_len));
        }
        match_len = 0;
      }
    }
    if (match_len > 0) {
      out.append(make_tuple(ms_1, ms_2, match_len));
    }
    return object(out);
  }

  double score;
  int offset;
  dna_sequence s1;
  dna_sequence s2;
  std::vector<align_state> m_vec;
};

namespace {

std::string dna_sequence_repr(const dna_sequence& seq) {
  return std::string("biograph.Sequence('") + seq.as_string() + std::string("')");
}

size_t dna_sequence_hash(dna_sequence seq) { return std::hash<std::string>{}(seq.as_string()); }

}  // namespace

void bind_dna_sequence(module& m) {
  class_<dna_sequence>(m, "Sequence",
                       "An object representing a DNA sequence. A sequence may be compared to "
                       "other\n"
                       "sequences and accessed in a manner similar to strings.\n"
                       "\n")
      .def(init<const std::string&>())
      .def(init<>())
      .def("__str__", &dna_sequence::as_string, "The nucleotide sequence as a simple string.\n")
      .def("__repr__", dna_sequence_repr,
           "The internal python representation of the sequence object.\n")
      .def("__len__", &dna_sequence::size, "The number of bases in the sequence.\n")
      .def("rev_comp", &dna_sequence::rev_comp,
           "Return a new sequence object of the reverse complement of "
           "this sequence.\n"
           "\n"
           "Example:\n"
           "\n"
           "    dna = Sequence('ACTG')\n"
           "    print dna.rev_comp() # CAGT\n")
      .def("__getitem__",
           [](const dna_sequence& seq, slice s) -> dna_sequence {
             ssize_t start, stop, step, slicelength;
             if (!(s.compute(seq.size(), &start, &stop, &step, &slicelength))) {
               throw error_already_set();
             }
             if (step != 1) {
               throw(std::out_of_range(
                   "Slicing a DNA sequence with steps other than 1 not supported"));
             }
             return seq.subseq(start, slicelength);
           })
      .def("__getitem__",
           [](const dna_sequence& seq, int pos) -> dna_sequence {
             if (pos < 0) {
               pos += int(seq.size());
             }
             if (pos >= int(seq.size()) || pos < 0) {
               throw(index_error());
             }
             return seq.subseq(pos, 1);
           })
      .def("__setitem__",
           [](dna_sequence& seq, const slice& s, const dna_sequence& replace_with) {
             ssize_t start, stop, step, slicelength;
             if (!(s.compute(seq.size(), &start, &stop, &step, &slicelength))) {
               throw error_already_set();
             }
             if (step != 1) {
               throw(index_error("Slicing a DNA sequence with steps other than 1 not supported"));
             }
             seq = seq.subseq(0, start) + replace_with + seq.subseq(stop, seq.size() - stop);
           })
      .def("__setitem__",
           [](dna_sequence& seq, const slice& s, const std::string& replace_with) {
             ssize_t start, stop, step, slicelength;
             if (!(s.compute(seq.size(), &start, &stop, &step, &slicelength))) {
               throw error_already_set();
             }
             if (step != 1) {
               throw(index_error("Slicing a DNA sequence with steps other than 1 not supported"));
             }
             seq = seq.subseq(0, start) + dna_sequence(replace_with) +
                   seq.subseq(stop, seq.size() - stop);
           })
      .def("__setitem__",
           [](dna_sequence& seq, int pos, char replace_with) {
             if (pos < 0) {
               pos += int(seq.size());
             }
             if (pos >= int(seq.size()) || pos < 0) {
               throw(index_error());
             }
             seq[pos] = dna_base(replace_with);
           })
      .def(self < self)
      .def(self <= self)
      .def(self == self)
      .def(self != self)
      .def(self > self)
      .def(self >= self)
      .def(self += self)
      .def(self + self)
      .def(
          "__radd__",
          [](const dna_sequence& seq, const std::string& str) -> dna_sequence { return str + seq; })
      .def("__hash__", dna_sequence_hash);

  implicitly_convertible<std::string, dna_sequence>();
}
