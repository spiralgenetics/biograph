"""
    core.py

    Core functest utility functions.
"""
import os
import os.path
import tempfile
import tools.py_version as py_version

# This tracks the current state of all executed tests. There must be a better way to do
# this with nose or unittest, but I can't find it.
TEST_RESULTS = dict()

for key in dir(py_version):
    if key.startswith('SPIRAL'):
        os.environ[key] = getattr(py_version, key).rstrip()

# Import from

search_path = os.environ["PATH"]
tmpdir = tempfile.mkdtemp()
os.mkdir(tmpdir + "/bin")
os.environ["PATH"] = tmpdir + "/bin:" + os.environ["PATH"]

def make_links(module_dir, base_binary, link_names):
    """Sets up links to module_directory/base_binary under the given link_names in $PATH"""
    if os.path.isfile(module_dir + "/" + base_binary):
        for alias in link_names:
            try:
                os.unlink(tmpdir + "/bin/" + alias)
            except OSError:
                pass
            os.symlink(os.getcwd() + "/" + module_dir + "/" + base_binary,
                       tmpdir + "/bin/" + alias)

make_links("modules/spec", "spec",
           ["spec", "fasta2ref", "ref2fasta", "szip", "bam2spec",
            "spec2bam", "spec_sample", "spec2fasta"])

make_links("modules/biograph", "bgbinary", ["bgbinary"])
