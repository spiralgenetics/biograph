#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_flat.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

// Generates a seqset for the given reads and their reverse
// complements.  Does not perform read correction.
std::unique_ptr<seqset_file> seqset_for_reads(const std::vector<dna_sequence>& reads);

// Generates a readmap for the given optionally-paired reads.  These
// reads must already exist in the seqset.  If present, saves the
// filename of the generated readmap in readmap_filename.
std::unique_ptr<readmap> readmap_for_reads(
    const std::shared_ptr<seqset>& the_seqset,
    const std::vector<std::pair<dna_sequence, dna_sequence>>& paired_reads,
    const std::vector<dna_sequence>& unpaired_reads,
    boost::optional<std::string*> readmap_filename = boost::none);

// Given a seqset, generates a seqset_flat containing the data
// therein.  Caller is responsible for making sure the_seqset outlives
// the produced seqset_flat.
std::unique_ptr<seqset_flat> seqset_flat_for_seqset(const seqset* the_seqset);

// Dumps a representation of the given seqset to stderr.  Useful for
// debugging.  Each line is prefixed by prefix.
void dump_seqset(const std::string& prefix, const seqset& the_seqset);

// Number of partitions to use when generating a test seqset.
extern unsigned g_seqset_build_partition_depth;

// Construct BioGraph seqset/readmap for vector of vector of dna_sequences.
// If unpaired, don't include the second dna_sequence in the vector, or make
// it's .size==0
// Paired/Unpaired can be provided together

std::pair<std::shared_ptr<seqset>, std::unique_ptr<readmap>> biograph_for_reads(
    const std::vector<std::vector<dna_sequence>>& all_reads);

MATCHER_P(SeqsetEntryIs, seq, "") {
  return ExplainMatchResult(seq, arg.sequence(), result_listener);
}

void PrintTo(const seqset_range& r, std::ostream* os);
