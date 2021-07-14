#include <pybind11/pybind11.h>
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "modules/io/utils.h"
#include "python/biograph/internal/assemble.h"
#include "python/biograph/internal/anchor.h"

variant_wrapper::variant_wrapper(const sv_out& var,
                                 std::shared_ptr<assembly> src,
                                 std::shared_ptr<reference> ref)
    : m_var(var),
      m_assembly(src),
      m_ref(ref),
      m_min_depth(255),
      m_max_depth(0),
      m_avg_depth(0.0) {
  for (int i = m_var.seq_begin - 1; i < m_var.seq_end + 1; i++) {
    uint8_t depth = m_assembly->depth[i];
    m_min_depth = std::min(depth, m_min_depth);
    m_max_depth = std::max(depth, m_max_depth);
    m_avg_depth += depth;
  }
  m_avg_depth /= 2 + m_var.seq_end - m_var.seq_begin;

  // canonicalize positions
  if (!left_forward() and !right_forward()) {
    m_var = m_var.flip(m_assembly->assembly.size());
  }
  // It can be reversed, just want left less than right
  if (left_position() > right_position()) {
    m_var = m_var.flip(m_assembly->assembly.size());
  }

}

std::string variant_wrapper::str() const {
  if (m_var.is_structural) {
    return printstring("%s:%lu%c_%d_%s:%lu%c", left_contig().c_str(),
                       left_position(), left_forward() ? '+' : '-',
                       m_var.seq_end - m_var.seq_begin, right_contig().c_str(),
                       right_position(), right_forward() ? '+' : '-');
  } else {
    return printstring(
        "%s:%lu_%lu:%d", left_contig().c_str(), left_position(),
        m_var.right_ref - m_var.left_ref, m_var.seq_end - m_var.seq_begin);
  }
}

std::string variant_wrapper::repr() const {
  return std::string("<biograph.Variant ") + str() + std::string(">");
}

std::string variant_wrapper::left_contig() const {
  size_t flat = m_var.left_ref.get_offset();
  return m_ref->get_assembly()
      .scaffold_order[m_ref->get_seq_position(flat).scaffold_id];
}

std::string variant_wrapper::right_contig() const {
  size_t flat = m_var.right_ref.get_offset();
  return m_ref->get_assembly()
      .scaffold_order[m_ref->get_seq_position(flat).scaffold_id];
}

unsigned long variant_wrapper::left_position() const {
  size_t flat = m_var.left_ref.get_offset();
  return m_ref->get_seq_position(flat).position;
}

unsigned long variant_wrapper::right_position() const {
  size_t flat = m_var.right_ref.get_offset();
  return m_ref->get_seq_position(flat).position;
}

dna_sequence variant_wrapper::variant_sequence() const {
  return m_assembly->assembly.subseq(m_var.seq_begin,
                                     m_var.seq_end - m_var.seq_begin);
}

dna_sequence variant_wrapper::assembly_sequence() const {
  return m_assembly->assembly;
}

reference_range variant_wrapper::range() const {
  size_t lflat = m_var.left_ref.get_offset();
  size_t rflat = m_var.right_ref.get_offset();
  return reference_range(m_ref, std::min(lflat, rflat),
                         std::max(lflat, rflat));
}

py_list variant_wrapper::depths() const {
  py_list ret;
  for (int i = m_var.seq_begin - 1; i < m_var.seq_end + 1; i++) {
    ret.append(m_assembly->depth[i]);
  }
  return ret;
}

py_list variant_wrapper::assembly_depths() const {
  py_list ret;
  for (uint8_t depth : m_assembly->depth) {
    ret.append(depth);
  }
  return ret;
}

variant_wrapper variant_wrapper::flip() const {
  variant_wrapper out;
  out.m_var = m_var.flip(m_assembly->assembly.size());
  out.m_assembly = std::make_shared<assembly>(m_assembly->flip());
  out.m_ref = m_ref;
  out.m_min_depth = m_min_depth;
  out.m_max_depth = m_max_depth;
  out.m_avg_depth = m_avg_depth;
  return out;
}

static void add_reference_depths(py_list out, dna_const_iterator a,
                                 dna_const_iterator b, const uint8_t* depths,
                                 const reference& ref) {
  bool reverse = false;
  if (a.is_rev_comp()) {
    dna_const_iterator t = a;
    a = b.rev_comp();
    b = t.rev_comp();
    reverse = true;
  }

  if (a.is_rev_comp() || b.is_rev_comp()) {
    throw io_exception("Impossible case #1 in add_reference_depths");
  }

  seq_position sp_a = ref.get_seq_position(a.get_offset());
  seq_position sp_b = ref.get_seq_position(b.get_offset());
  if (sp_a.scaffold_id != sp_b.scaffold_id) {
    throw io_exception("Impossible case #2 in add_reference_depths");
  }

  py_list py_dp;
  for (size_t i = 0; i < sp_b.position - sp_a.position + 1; i++) {
    if (reverse) {
      py_dp.append(depths[sp_b.position - sp_a.position - i]);
    } else {
      py_dp.append(depths[i]);
    }
  }

  std::string scaf = ref.get_assembly().scaffold_order[sp_a.scaffold_id];
  // printf("REF: %s: [%lu, %lu]\n", scaf.c_str(), sp_a.position,
  // sp_b.position);

  py_list all;
  all.append(scaf);
  all.append(sp_a.position + 1);
  all.append(py_dp);
  out.append(all);
}

py_list py_assemble(py_list src, py_list dest, uint8_t min_overlap,
                    uint32_t max_steps, bool skip_ambig,
                    const std::shared_ptr<readmap>& wrap) {
  std::shared_ptr<seqset_bitmap_base> bm = wrap;

  py_list ret;
  py_list ref_seq;

  if (pybind11::len(src) == 0 || pybind11::len(dest) == 0) {
    return ret;
  }

  anchor_wrapper x = pybind11::cast<anchor_wrapper>(src[0]);
  anchor_wrapper y = pybind11::cast<anchor_wrapper>(dest[0]);

  if (x.m_reference != y.m_reference || x.m_seqset != y.m_seqset) {
    throw io_exception("Using anchors with mismatched References or Seqsets");
  }

  const reference& ref = *x.m_reference;

  auto out_it =
      boost::make_function_output_iterator([&](const assembly& assem) {

        std::shared_ptr<assembly> pa = std::make_shared<assembly>(assem);
        dna_slice left_range = ref.get_supercontig(assem.left.get_offset());
        dna_slice right_range = ref.get_supercontig(assem.right.get_offset());

        if (assem.left.is_rev_comp()) {
          left_range = left_range.rev_comp();
        }

        if (assem.right.is_rev_comp()) {
          right_range = right_range.rev_comp();
        }

        left_range = dna_slice(assem.left, left_range.end());
        right_range = dna_slice(right_range.begin(), assem.right + 1);

        std::vector<sv_out> bits =
            call_structural(assem.assembly, left_range, right_range, 200);
        if (bits.size() == 0) {
          return;
        }

        std::shared_ptr<assembly> par;
        dna_const_iterator it_ref = assem.left;
        const uint8_t* depth_start = &assem.depth[0];

        for (size_t i = 0; i < bits.size(); i++) {
          if (bits[i].anchor_drop) {
            // printf("Drop in assembly: %d\n", assem.id);
            if (bits[i].seq_begin == 0) {
              it_ref = bits[i].right_ref;
              pa->left = it_ref;
            } else {
              pa->right = bits[i].left_ref;
            }
            continue;
          }

          add_reference_depths(ref_seq, it_ref, bits[i].left_ref, depth_start,
                               ref);
          it_ref = bits[i].right_ref;
          depth_start = &assem.depth[bits[i].seq_end];

          // Get Ref between it_ref + left_ref
          // Update it_ref to right_ref
          if (!bits[i].is_structural && bits[i].left_ref.is_rev_comp()) {
            // TODO: This code is UNTESTED!!!
            sv_out rb = bits[i];
            rb.left_ref = bits[i].right_ref.rev_comp();
            rb.right_ref = bits[i].left_ref.rev_comp();
            rb.seq_begin = assem.assembly.size() - bits[i].seq_end;
            rb.seq_end = assem.assembly.size() - bits[i].seq_begin;
            if (!par) {
              par = std::make_shared<assembly>(pa->flip());
            }
            ret.append(variant_wrapper(rb, par, x.m_reference));
          } else {
            ret.append(variant_wrapper(bits[i], pa, x.m_reference));
          }
        }

        add_reference_depths(ref_seq, it_ref, pa->right, depth_start, ref);
      });

  auto transform_func = [&](pybind11::handle h) {
    return pybind11::cast<const anchor_wrapper&>(h).m_anchor;
  };

  auto left = boost::make_iterator_range(
      boost::make_transform_iterator(src.begin(), transform_func),
      boost::make_transform_iterator(src.end(), transform_func));

  auto right = boost::make_iterator_range(
      boost::make_transform_iterator(dest.begin(), transform_func),
      boost::make_transform_iterator(dest.end(), transform_func));

  seqset_assemble(out_it, x.m_seqset->get_seqset(), left, right, min_overlap,
                  max_steps, skip_ambig, *wrap);

  py_list both;
  both.append(ret);
  both.append(ref_seq);
  return both;
}

bool variant_wrapper::py_lessthan(variant_wrapper other) {
  unsigned long mStart =
      std::min(this->left_position(), this->right_position());
  unsigned long oStart =
      std::min(other.left_position(), other.right_position());
  return mStart < oStart;
}

int variant_wrapper::py_cmp(variant_wrapper other) {
  if (this->left_contig() == other.left_contig() &&
      this->right_contig() == other.right_contig() &&
      this->left_position() == other.left_position() &&
      this->right_position() == other.right_position() &&
      this->variant_sequence() == other.variant_sequence()) {
    return 0;
  }
  if (this->py_lessthan(other)) {
    return -1;
  }

  return 1;
}

using namespace pybind11;
void bind_assemble(module& m) {
  class_<variant_wrapper>(m,
      "Variant",
      "An object representing a variation from reference.\n"
      "\n"
      "Left and Right correspond to where assemblies that contain variation from the reference align\n"
      "For example, a deletion's assembly would look like the following:\n"
      "  ATCAG----CAGAT\n"
      "  ATCAGAGATCAGAT\n"
      "Where:\n"
      "  variant.left_position = 4\n"
      "  variant.right_position = 9\n"
      "  variant.sequence = ''\n"
      "  variant.range.sequence = 'GAGAT'\n"
      "For indels, left_forward and right_forward will not match (one True,\n"
      "one False).\n"
      "\n"
      "Variants with both left_forward and right_forward True or False are\n"
      "inversions.\n"
      "")
      .def("__str__", &variant_wrapper::str,
           "A unique identifier for this variant.")
      .def("__repr__", &variant_wrapper::repr,
           "The python representation of this variant object.\n")
      .def("__lt__", &variant_wrapper::py_lessthan,
           "Is this variant's chrom:start less than another's\n")
      .def("__cmp__", &variant_wrapper::py_cmp,
           "Compare these variants' position and sequence")
      .def_property_readonly("is_structural", &variant_wrapper::is_structural,
                    "True if this is a structural variant.")
      .def_property_readonly(
           "left_contig", &variant_wrapper::left_contig,
           "The left side contig for this variant. Structural variants may\n"
           "cross supercontig boundaries.")
      .def_property_readonly(
           "right_contig", &variant_wrapper::right_contig,
           "The right side contig for this variant. Structural variants may\n"
           "cross supercontig boundaries.")
      .def_property_readonly("left_position", &variant_wrapper::left_position,
                    "The left position for this variant (pointing to the "
                    "reference base\n"
                    "immediately before the variant).")
      .def_property_readonly("right_position", &variant_wrapper::right_position,
                    "The right position for this variant (pointing to the "
                    "reference base\n"
                    "immediately after the variant).")
      .def_property_readonly("left_forward", &variant_wrapper::left_forward,
                    "True if the left side is in the forward direction. Always "
                    "true for\n"
                    "non-structural variants.")
      .def_property_readonly("right_forward", &variant_wrapper::right_forward,
                    "True if the right side is in the forward direction. "
                    "Always true for\n"
                    "non-structural variants.")
      .def_property_readonly("sequence", &variant_wrapper::variant_sequence,
                    "A Sequence object representing the variant nucleotide\n"
                    "sequence (if any).")
      .def_property_readonly("range", &variant_wrapper::range,
                    "A ReferenceRange over the [self.left_position:self.right_position) in\n"
                    "the reference where this variant occurs.")
      .def_property_readonly(
           "depths", &variant_wrapper::depths,
           "A list of read depth counts for each base in this variant,\n"
           "including the bases immediately preceding and following it in\n"
           "the reference.\n"
           "\n"
           "For example, a SNP would have depth list containing three values:\n"
           "the count of reads that match the last base in the reference "
           "prior\n"
           "to the SNP, the SNP itself, and the next base in the reference.")
      .def_property_readonly(
           "min_depth", &variant_wrapper::min_depth,
           "The minimum depth of coverage of any base in this variant.")
      .def_property_readonly(
           "max_depth", &variant_wrapper::max_depth,
           "The maximum depth of coverage of any base in this variant.")
      .def_property_readonly(
           "avg_depth", &variant_wrapper::avg_depth,
           "The average depth of coverage of all bases in this variant (float)")
      .def_property_readonly("min_overlap", &variant_wrapper::min_overlap,
                    "The lowest overlap found in the entire assembly.")
      .def_property_readonly(
           "assembly_id", &variant_wrapper::assembly_id,
           "A unique identifier for the assembly where this variant was "
           "found.\n"
           "Note that several variants may be found within a single assembly.")
      .def_property_readonly("assembly_sequence", &variant_wrapper::assembly_sequence,
                    "A Sequence object representing the entire assembly\n"
                    "including all variants.")
      .def_property_readonly("assembly_depths", &variant_wrapper::assembly_depths,
                    "A list of coverage depths for every base in the assembled "
                    "sequence.")
      .def_property_readonly("assembly_begin", &variant_wrapper::assembly_begin,
                    "The offset within this assembly where the variant begins "
                    "(inclusive).")
      .def_property_readonly("assembly_end", &variant_wrapper::assembly_end,
                    "The offset within this assembly where the variant ends "
                    "(exclusive).")
      .def("flip", &variant_wrapper::flip,
           "A Sequence object for the reverse complement of this variant.");

  m.def("assemble", py_assemble,
      "assemble([src], [dest], min_overlap, max_steps, skip_ambiguous, "
      "readmap)\n"
      "\n"
      "Find all possible assemblies between the lists of anchors in src\n"
      "and dest, for all reads matching the given Readmap.\n"
      "\n"
      "Stop when all assemblies are found or max_steps have\n"
      "been taken. If skip_ambiguous is True, don't assemble across\n"
      "ambiguous regions of the reference.\n"
      "\n"
      "Returns a list containing two sublists: \n"
      "\n"
      "  * assemblies[0] contains the variants called from all assemblies\n"
      "    found, with one entry for each variant. Since there can be\n"
      "    many variants for each assembly, this list can be much larger\n"
      "    than the total number of assemblies.\n"
      "\n"
      "  * assemblies[1] contains a list of each assembly's coverage\n"
      "    depth. There is one entry for each assembly found.\n"
      "\n"
      "NOTE: Some of the coverage of each assembly contributes to the\n"
      "reference coverage (ie. the parts of the anchor that match\n"
      "reference). This coverage is not included in reference_range\n"
      "coverage, since that can only include reads that precisely\n"
      "match reference.\n"
      "\n"
      "To calculate the true reference coverage for the anchors,\n"
      "you will need to add the assembly coverage for bases that\n"
      "lie in the reference back to the ReferenceRange coverage.\n"
      "\n"
      "Each coverage list is of the format:\n"
      "\n"
      "    ['supercontig', position, [25,25,26,26,26...]]\n"
      "\n"
      "Example:\n"
      "\n"
      "  from biograph import BioGraph, Reference, find_anchors, assemble\n"
      "  my_bg = BioGraph('NA12878.bg')\n"
      "  ref = Reference('/reference/human_g1k_v37/')\n"
      "  rr = ref.make_range('2', 1, 200000, False)\n"
      "  fwd = find_anchors(my_bg, rr, True, 70, 10000)\n"
      "  rev = find_anchors(my_bg, rr, False, 70, 10000)\n"
      "  both = assemble(fwd, rev, 70, 100000, True)\n"
      "  print both[0][0].sequence # the first variant\n");
}
