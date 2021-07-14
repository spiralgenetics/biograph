#pragma once

#include <string>
#include <memory>

#include "modules/bio_base/flat_ref.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/dna_sequence.h"

// Creates a flat ref with the given sequences.  Each will be in its
// own scaffold, and will start at offset 0 within that scaffold.
std::unique_ptr<flat_ref> create_flat_ref(std::vector<dna_sequence> extents);

// create_reference returns a full reference complete with BWT
// containing the given extents.  It also adds extents of length 1 for
// each base to avoid errors building the BWT when at least 1 of each
// base isn't present.
std::unique_ptr<reference> create_reference(const std::vector<dna_sequence>& extents);

// This version of create_reference allows "N"s to be included in the
// reference.  Note that small sequences of "N"s may be silently
// converted to arbitrary bases, instead of separating extents.  See
// fasta_ref_importer::add_base.
std::unique_ptr<reference> create_reference_str(const std::vector<std::string>& extents);

// Allow gtest to pretty print values of type 'flat_ref::extent_t'
inline void PrintTo(const flat_ref::extent_t& extent, ::std::ostream* os) {
  (*os) << "\n  Scaffold: " << extent.scaffold_name << " Offset: " << extent.offset
        << " Size: " << extent.size << " Flat: " << extent.flat;
}

// Allow gtest to pretty print values of type 'spec_header::scaffold_t':
inline void PrintTo(const spec_header::scaffold_t& scaffold, ::std::ostream* os) {
  (*os) << "\n  Name: " << scaffold.name << " size: " << scaffold.size;
}
