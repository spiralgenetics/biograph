"""Loads a dynamic module by the name of MODULE_NAME_ext.  Used from a
spiral_py_dynamic_module rule.
"""

# pylint:disable=undefined-variable,bare-except

import imp
import os
import os.path
import sys
import tempfile

def _get_dyn_module_contents():
    """Returns the contents of the dynamic module .so file as a string"""
    dirname = os.path.dirname(__file__)
    filename = dirname + "/MODULE_NAME_internal.so"

    # First, try to load from the loader that loaded this .py file, in
    # case it's e.g. a zipimporter:
    try:
        return __loader__.get_data(filename)
    except:
        pass

    # Otherwise, try to load the file from the filesystem.
    try:
        with open(filename) as f:
            return f.read()
    except:
        pass

    raise ImportError("Unable to find dynamic library " + filename)


def _import_dynamic_module():
    """Imports the dynamic module, and replaces this module with it."""
    file_contents = _get_dyn_module_contents()

    (fdnum, temp_filename) = tempfile.mkstemp(suffix=".so")
    f = os.fdopen(fdnum, "w")
    f.write(file_contents)
    f.close()

    # save unlink since we're replacing the module we're currently
    # executing in.
    unlink = os.unlink
    try:
        # Replace current module with dynamically loaded one
        new_module = imp.load_dynamic("MODULE_NAME", temp_filename)
        sys.modules[__name__] = new_module
    finally:
        unlink(temp_filename)

_import_dynamic_module()
