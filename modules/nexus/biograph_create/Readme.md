<!-- dx-header -->
# BioGraph: Create/Variants (DNAnexus Platform App)

Create a new BioGraph from reads and call variants against a reference.
<!-- /dx-header -->

# Software License

This is part of the BioGraph v4.0.0 limited license.

# What does this app do?

This app will run the BioGraph Engine and index/compress your input BAM/CRAM file into the BioGraph Format by running `biograph create`. The resulting BioGraph will then be compared against your input Reference file using the `biograph variants` command to create a VCF output of all discovered alleles.

# What inputs are required?

* `input_reads` - The input set of sequencing reads to be processed in BAM or CRAM format
* `reference` - A reference genome archive created by the `BioGraph: reference` app
* `sample_id` - A unique name to identify this sample

# What does this app output?

The app will create a top level directory in the current project, named for the sample_id provided. For example, if the sample_id is NA12878, the following files will be created:

* `NA12878.bg/seqset` - The compressed/indexed minimal set of maximal kmers
* `NA12878.bg/variants.vcf.gz` - Sorted/compressed variant call format file of every allele discovered by `biograph variants`
* `NA12878.bg/variants.vcf.gz.tbi` - Tabix index of the `*.vcf.gz`
* `NA12878.bg/coverage/*.readmap` - Read coverage data
* `NA12878.bg/metadata/bg_info.json` - Json file containing metadata about this BioGraph
* `NA12878.bg/qc/create_log.txt` - Run details for the `create` command
* `NA12878.bg/qc/create_stats.json` - Key/Value pairs of QC stats collected during the BioGraph create stage
* `NA12878.bg/qc/kmer_quality_report.html` - Interactive report of kmerization and error correction stats generated during `biograph create`
* `NA12878.bg/qc/variants_log.txt` - Run details for the `variants` command
* `NA12878.bg/qc/variants_stats.json` - Key/Value pairs of QC stats collected during the BioGraph variants stage

# How does this app work?

See the `biograph create` and `biograph variants` commands in the BioGraph User Guide at https://www.spiralgenetics.com/user-documentation for details.

Additional options can be input in the "Advanced 'biograph create' options" and "Advanced 'biograph variants' options" found in the App.
