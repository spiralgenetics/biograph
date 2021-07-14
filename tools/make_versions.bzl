load(
    "//:versions.bzl",
    "SPEC_VERSION",
    "BIOGRAPH_VERSION",
    "SEQSET_VERSION",
)

VERSION_TABLE = {
    "SPEC_VERSION": SPEC_VERSION,
    "BIOGRAPH_VERSION": BIOGRAPH_VERSION,
    "SEQSET_VERSION": SEQSET_VERSION,
}

def make_version_h(name, out_file):
    """Generates a .h file containing #defines for all the versions in VERSION_TABLE."""
    cmd = """(
echo '#pragma once'
"""
    for version_id, version in VERSION_TABLE.items():
        cmd = cmd + "echo '#define " + version_id + " \"" + version + "\"'\n"
    cmd = cmd + """
) > $@"""
    native.genrule(name=name, outs=[out_file],
                   cmd=cmd)

def make_version_py(name, out_file):
    """Generates a .py file containing version definitions for all the
    versions in VERSION_TABLE."""

    cmd = "("
    for version_id, version in VERSION_TABLE.items():
        cmd = cmd + "echo '" + version_id + "=\"" + version + "\"'\n"
    cmd = cmd + ") > $@"

    native.genrule(name=name, outs=[out_file], cmd=cmd)
