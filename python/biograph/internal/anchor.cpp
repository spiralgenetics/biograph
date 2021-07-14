
#include "python/biograph/internal/anchor.h"
#include <boost/function_output_iterator.hpp>

anchor_wrapper::anchor_wrapper(const anchor& anchor, std::shared_ptr<seqset_file> seqset,
                               std::shared_ptr<reference> reference)
    : m_anchor(anchor), m_seqset(seqset), m_reference(reference) {}

seqset_range anchor_wrapper::read() const {
  return m_seqset->get_seqset().ctx_entry(m_anchor.entry);
}

reference_range anchor_wrapper::range() const {
  SPLOG("Overlap = %d", m_anchor.overlap);
  size_t flat_end = m_anchor.ref_pos.get_offset();
  size_t flat_start = (m_anchor.ref_pos + (m_anchor.overlap - 1)).get_offset();
  return reference_range(m_reference, std::min(flat_start, flat_end),
                         std::max(flat_start, flat_end) + 1);
}

bool anchor_wrapper::forward() const { return m_anchor.ref_pos.is_rev_comp(); }

// TODO: make forward an enum for direction (FORWARD, REVERSE, BOTH)
py_list py_anchor(const std::shared_ptr<biograph>& my_bg, const reference_range& ref_range,
                  bool forward, uint8_t min_overlap, uint32_t max_anchors) {
  py_list r;
  std::shared_ptr<seqset_file> seqset = my_bg->get_seqset();
  std::shared_ptr<readmap> rm = my_bg->open_readmap("");
  std::shared_ptr<seqset_bitmap_base> bm = rm;
  std::shared_ptr<reference> ref = ref_range.get_reference();
  dna_slice ref_slice(ref->get_dna(ref_range.get_flat_start()),
                      ref->get_dna(ref_range.get_flat_end()));
  if (forward) {
    ref_slice = ref_slice.rev_comp();
  }
  seqset_anchor(boost::make_function_output_iterator(
                    [&](const anchor& a) { r.append(anchor_wrapper(a, seqset, ref)); }),
                seqset->get_seqset(), ref_slice, min_overlap, max_anchors, *bm);
  return r;
}

using namespace pybind11;
void bind_anchor(module& m) {
  class_<anchor_wrapper>(m,"Anchor",
                         "An Anchor object is a read that partially matches the reference.\n"
                         "")
      .def_property_readonly("read", &anchor_wrapper::read,
                    "A SeqsetEntry pointing to the read that contains this anchor.\n")
      .def_property_readonly("range", &anchor_wrapper::range,
                    "A ReferenceRange object for the region of reference that overlaps\n"
                    "this anchor.")
      .def_property_readonly("forward", &anchor_wrapper::forward,
                    "True if the anchor departs from reference in the forward (5'->3')\n"
                    "direction, otherwise False.\n"
                    "\n"
                    "Trivial example, assuming min_overlap of 3:\n"
                    "\n"
                    "               ACTATGATGC  << anchor in the forward direction\n"
                    "                 TATGATGC  << not returned, less specific\n"
                    "  ref: ...CTTGAACTATGATG...\n"
                    "         GGAACTTG          << anchor in the reverse direction\n")
      .def_property_readonly("overlap", &anchor_wrapper::overlap,
                    "The number of bases in this Anchor that match the reference.\n");

  m.def("find_anchors", py_anchor, arg("biograph"), arg("ref_range"), arg("forward") = true,
        arg("min_overlap") = 70, arg("max_anchors") = 10000,
        "find_anchors(biograph, ref_range, forward=True, min_overlap=70, "
        "max_anchors=10000)\n"
        "\n"
        "Return a list of Anchors for a BioGraph in the given Reference "
        "Range.\n"
        "Context is ignored.\n"
        "\n"
        "The forward parameter is a Boolean to specify whether to look for\n"
        "Anchors in the forward (True) or reverse (False) direction. "
        "Anchoring\n"
        "stops at the end of the ref_range or when max_anchors have been "
        "found.\n"
        "\n"
        "This function returns a list of reads that have more than\n"
        "min_overlap bases in common with reference in that range. Note\n"
        "that Anchor() only returns the MOST specific anchors for the\n"
        "given query; otherwise valid Anchors that are a subset of a larger\n"
        "Anchor are not returned.\n"
        "\n"
        "Example:\n"
        "\n"
        "    from biograph import BioGraph, Reference, find_anchors\n"
        "    my_bg = BioGraph('NA12878_S1.seqset', 'NA12878_S1.readmap')\n"
        "    ref = Reference('/reference/human_g1k_v37/')\n"
        "    rr = ref.make_range('2', 1, 200000, False)\n"
        "    fwd = find_anchors(my_bg, rr, True, 70, 10000)\n"
        "    rev = find_anchors(my_bg, rr, False, 70, 10000)\n");
}
