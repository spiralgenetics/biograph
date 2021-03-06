<!-- dx-header -->
# BioGraph (DNAnexus Platform App)

Spiral's BioGraph full_pipeline v6.0.1

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.
<!-- /dx-header -->

<!-- Insert a description of your app here -->
This app is a wrapper around the `biograph full_pipeline` command.
It will take a single-sample's bam or cram file, and produce the BioGraph file and the results.vcf.gz

This is different from documented BioGraph's pipeline in that:

1. There is no access to `--resume` or `--stop` for subsetting the pipeline.
2. The produced BioGraph is packaged into a `.tar` file for ouptut
3. The results.vcf.gz and its index are not kept inside the BioGraph file's `analysis/` directory.
4. The input BioGraph reference is a `.tar` file that will be un-tarred after download

!!IMPORTANT!!

It is HIGHLY recommended that analysis is limited to the main chromosomes for Human genomes (e.g. chr1-22, X, Y)
To do this restriction, set the `--discovery` parameter's value to `--bed reference/grch38/regions.bed`, where "grch38"
is the name of the reference used.

<!--
TODO: This app directory was automatically generated by dx-app-wizard;
please edit this Readme.md file to include essential documentation about
your app that would be helpful to users. (Also see the
Readme.developer.md.) Once you're done, you can remove these TODO
comments.

For more info, see https://documentation.dnanexus.com/developer.
-->
