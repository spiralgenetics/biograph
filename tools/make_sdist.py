#!/usr/bin/python3

"""Creates a "sdist" setuptools package for the BioGraph SDK.

Many of the parameters of this script are hardcoded to support the
Biograph SDK.  If this script needs to package other things, it
should be pretty striaghtforward to parameterize package names and
descriptions.
"""

from __future__ import print_function
import tempfile
import argparse
import importlib
import os
import os.path
import subprocess
import re

import sys
import shutil
from glob import glob

PACKAGE_FILE_REGEX = re.compile(".*\\.(py|so)$")

parser = argparse.ArgumentParser(description="Process some integers.")
parser.add_argument("packages", nargs="+", help="Packages to include")
parser.add_argument("--script", dest="scripts", action="append", default=[],
                    help="Path of a script to install in 'bin'; may be specified multiple times")
parser.add_argument(
    "--readme", help="Path to a 'README' file to use; name an extension are preserved")
parser.add_argument("--out", dest="output_file", required=True,
                    help="Output filename to store the generated .tar.gz package in")
parser.add_argument("--bin", dest="bins", action="append", default=[],
                    help="Compilde binary executables that should be included in the packgae.")
parser.add_argument("--requirements", dest="requirements",
                    help="Requirements file")
parser.add_argument(
    "--verbose", "-v", help="Verbose output.", action="store_true")


def pretty_repr(source):
    """Formats similar to repr() but takes advatage of triple quotes."""
    if not "\n" in source:
        # Short enough that triple quotes aren't particularly useful.
        return repr(source)
    for triple_quotes in ['"""', "'''"]:
        if triple_quotes not in source:
            return triple_quotes + source + triple_quotes
    # Unable to use triple quotes, since it already contains both kinds.
    return repr(source)

args = parser.parse_args()
dist_dir = tempfile.mkdtemp()
if args.verbose:
    print("Using dist dir %s" % dist_dir)

version = None
for pkgname in args.packages:
    if args.verbose:
        print("Processing package " + pkgname)
    dirname = pkgname.replace(".", "/")
    if not os.path.isdir(dist_dir + "/" + dirname):
        os.makedirs(dist_dir + "/" + dirname)
    pkgobj = importlib.import_module(pkgname)
    if not version:
        try:
            version = pkgobj.version()
        except AttributeError:
            pass

    seen_files = dict()

    for path in pkgobj.__path__:
        if args.verbose:
            print("Processing path " + path)
        for filename in os.listdir(path):
            if not re.match(PACKAGE_FILE_REGEX, filename):
                continue
            if filename in seen_files:
                print("Already saw {filename} in {old_path}, but also present in {new_path}".
                      format(filename=filename, old_path=seen_files[filename], new_path=path))
                sys.exit(1)
            seen_files[filename] = path
            if args.verbose:
                print("About to symlink " + filename)
            os.symlink(
                path + "/" + filename, dist_dir + "/" + dirname + "/" + filename)

if args.scripts or args.bins:
    os.makedirs(dist_dir + "/bin")
    for script in (args.scripts + args.bins):
        os.symlink(os.path.abspath(script),
                   dist_dir + "/bin/" + os.path.basename(script))

long_description = ""
if args.readme:
    os.symlink(os.path.abspath(args.readme),
               dist_dir + "/" + os.path.basename(args.readme))
    long_description = open(args.readme).read()

if not version:
    print("Unable to determine SDK version")
    sys.exit(1)

reqs = []
if args.requirements:
    with open(args.requirements) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            reqs.append(line)

package_list = ",".join(['"' + pkg + '"' for pkg in args.packages])
script_list = ",".join(
    ['"bin/' + os.path.basename(script) + '"' for script in args.scripts])
bin_list = ",".join(
    ['("bin", ["bin/' + os.path.basename(binfile) + '"])' for binfile in args.bins])
requires_list = ",".join('"' + req + '"' for req in reqs)
with open(dist_dir + "/setup.py", "w") as setup_file:
    setup_file.write("""
#!/usr/bin/python3
'''For distributing and installing BioGraph'''

from setuptools import setup

setup(
    name='BioGraph',
    version='""" + version + """',
    packages=[""" + package_list + """],
    author='SpiralGenetics',
    author_email='support@spiralgenetics.com',
    scripts=[""" + script_list + """],
    license='LICENSE.txt',
    url="www.spiralgenetics.com",
    description='BioGraph genomics pipeline',
    long_description=""" + pretty_repr(long_description) + """,
    package_data={'': ["*.so"]},
    data_files=[""" + bin_list + """],
    zip_safe=False,
    entry_points={
      'console_scripts': [
         'biograph = biograph.__main__:main'
      ]
    },
    install_requires=[""" + requires_list + """]
)
""")

subprocess.check_call(
    "cd " + dist_dir + " && python3 setup.py -q sdist", shell=True)
matched_files = glob(dist_dir + "/dist/*.tar.gz")
if len(matched_files) != 1:
    print("%d unexpected dist files generated: %s" %
          (len(matched_files), repr(matched_files)))
    sys.exit(1)
shutil.move(matched_files[0], args.output_file)
