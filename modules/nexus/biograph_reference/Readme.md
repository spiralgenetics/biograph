<!-- dx-header -->
# BioGraph: Reference (DNAnexus Platform App)

Create a new BioGraph Reference archive from a FASTA file.
<!-- /dx-header -->

# Software License

This is part of the BioGraph v4.0.0 limited license.

# What does this app do?

This app will convert a FASTA file to the BioGraph Reference archive format, for use with the BioGraph: create app.

# What inputs are required?

* `fasta` - A genomic reference in FASTA format, optionally gzipped
* `reference_name` - A name for this reference (eg. hs37d5)

# What does this app output?

* `reference/*.tgz` - A gzipped tar archive of the BioGraph reference directory

# How does this app work?

See the `biograph reference` command in the BioGraph User Guide at https://www.spiralgenetics.com/user-documentation
