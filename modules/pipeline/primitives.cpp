#include "tools/version.h"
#include "modules/pipeline/primitives.h"
#include "modules/pipeline/steptype.h"
#include "modules/pipeline/all_steps.h"

void add_primitives()
{
	// Data types

	datatype_registry::rest_register("/datatypes");

	datatype_registry::add(std::make_shared<datatype>("unaligned_reads",
		"Unaligned Reads",
		"Sequence read basecalls and quality values"));
	datatype_registry::add(std::make_shared<datatype>("qseq", "qseq", "qseq"));
	datatype_registry::add(std::make_shared<datatype>("fasta", "fasta", "fasta"));
	datatype_registry::add(std::make_shared<datatype>("unaligned_read_set",
		"Unaligned Reads (deprecated)",
		"Sequence read basecalls and quality values (deprecated)"));
	datatype_registry::add(std::make_shared<datatype>("bam",
		"BAM",
		"Binary SAM files"));
	datatype_registry::add(std::make_shared<datatype>("sequence",
		"Sequence",
		"DNA sequence data"));
	datatype_registry::add(std::make_shared<datatype>("read_qual",
		"Read Quality Report",
		"Quality statistics on a set of unaligned reads"));
	datatype_registry::add(std::make_shared<datatype>("reference",
		"Reference Database",
		"Reference genome data"));
	datatype_registry::add(std::make_shared<datatype>("variants",
		"Variation Data",
		"Individual variation data"));
	datatype_registry::add(std::make_shared<datatype>("kmer_count",
		"Kmer Counts",
		"Kmers and their abundance"));
	datatype_registry::add(std::make_shared<datatype>("kmer_scores",
		"Kmer Scores",
		"Kmers and kmer quality scores"));
	datatype_registry::add(std::make_shared<datatype>("fusion_count",
		"Fusion Counts",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("covar_table",
		"Base Quality Recalibration Covariate Table",
		"Table generated by first step of base quality recalibration"));
	datatype_registry::add(std::make_shared<datatype>("kmer_db",
		"Kmer Database",
		"Database generated by kmer counts"));
	datatype_registry::add(std::make_shared<datatype>("kmer_aligned_reads",
		"Kmer Aligned Reads",
		"Reads ready to be called into variants"));
	datatype_registry::add(std::make_shared<datatype>("recal_qual",
		"Base Quality Recalibration",
		"Second step of base quality recalibration that adjusts the alignment qualities"));
	datatype_registry::add(std::make_shared<datatype>("realign_indels",
		"Indel local realignment",
		"Locally realign near indels"));
	datatype_registry::add(std::make_shared<datatype>("corrected_reads",
		"Corrected reads",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("anchored_reads",
		"Anchored reads",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("histogram",
		"Histogram of some data",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("struct_var",
		"Structual variation data",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("read_support",
		"Structural variation read support data",
		"Reads used to generate structural variation data"));
	datatype_registry::add(std::make_shared<datatype>("pbwt",
		"Prefix BWT table",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("coverage_records",
		"Coverage records",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("coverage_db",
		"Coverage database",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("sv_calls",
		"Called sv data",
		"Some description here"));
	datatype_registry::add(std::make_shared<datatype>("prefixes",
		"Prefixes",
		"Some description here"));

	// Steps

	steptype_registry::rest_register("/steptypes");
}
