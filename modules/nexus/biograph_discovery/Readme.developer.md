# BioGraph Developer Readme

<!--
TODO: Please edit this Readme.developer.md file to include information
for developers or advanced users, for example:

* Information about app internals and implementation details
* How to report bugs or contribute to development
-->

To build this app:
1. Use `dx build biograph_discovery` while logged into you your account and the
appropriate project.
2. You can run either from the DNANexus web interface or via:

```
dx run biograph_full_pipeline \
	-itarball=BioGraph-6.0.4.tgz \
	-ibiograph=biograph_file.tgz \
	-ireference=hg19.tar \
	-imodel=biograph_model.ml \
	-ireads=some.cram \
	-idiscovery="--bed reference/hg19/regions.bed" \
	--destination=output_dir/
```

IMPORTANT!!

1. The BioGraph tarball is uploaded as input into the job.
2. The BioGraph reference need to be tarballs.
3. The BioGraph reference is extracted into a directory named `reference/`
4. The created BioGraph Format file will be saved as a tarball.
5. The resultant VCF will be saved beside the created BioGraph, not inside. There's a parameter to switch that.
6. There's a parameter to not upload the BioGraph Format tarball that's useful for `--resume` operations.

## Running this app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See <a
href="https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification">Run
Specification</a> in the API documentation for more information about the
available instance types.
