"""Runs 'doctest' examples in all biograph modules and submodules"""

import sys
import doctest

import biograph
import biograph.variants

# Boost sets the __module__ attribute of class methods to be the name
# class they're contained in, which is wrong.  Instead, it should
# store this in __objclass__.  But anyways, that means that doctest
# thinks that class methods are part of a different module and doesn't
# search them.

VERBOSE = False
DISPLAY_PRUNED_OBJECTS = False

MODULES_TO_TEST = (
    biograph,
    biograph.utils,
    biograph.variants,
    biograph._capi, # pylint:disable=protected-access
)

boost_builtin_function_or_method = biograph._capi.version.__class__ # pylint:disable=protected-access
did_prune = dict()
test_modules_names = {mod.__name__ for mod in MODULES_TO_TEST}

def patch_from_module(orig_from_module):
    """Monkey-patch the given doctest.DocTestFinder._from_module"""
    # TODO(nils): This is a kludge; we should either fix boost or find some other
    # way to solve this problem.
    def new_from_module(doctest_self, module, obj):
        """Monkey-patched replacement for doctest.DocTestFinder._from_module"""
        # If it's boost's "builtin_function_or_method", always
        # recurse into it no matter what its __module__ says.
        if isinstance(obj, boost_builtin_function_or_method):
            return True
        result = orig_from_module(doctest_self, module, obj)
        if result or not DISPLAY_PRUNED_OBJECTS:
            return result
        if id(obj) not in did_prune:
            obj_name = "(unknown)"
            obj_module = "(unknown)"
            try:
                obj_name = obj.__name__
            except AttributeError:
                pass
            try:
                obj_module = obj.__module__
            except AttributeError:
                pass
            if obj_module not in test_modules_names:
                print("Skipping out-of-module {} which says its name is {} and is in {}".format(
                    obj, obj_name, obj_module))
            did_prune[obj] = True
        return False
    return new_from_module

doctest.DocTestFinder._from_module = patch_from_module(doctest.DocTestFinder._from_module) # pylint:disable=protected-access

globs = {
    "my_ref": biograph.Reference('datasets/lambdaToyData/benchmark/ref_lambda'),
    "human_ref": biograph.Reference('/reference/human_g1k_v37'),
    "my_biograph": biograph.BioGraph('datasets/lambdaToyData/benchmark/family_lambda.bg'),
    "biograph" : biograph,
}

failed_modules = []

def test_mod(mod):
    """Runs doctest tests on the given module"""
    modname = mod.__name__
    fail, total = doctest.testmod(mod, verbose=VERBOSE, globs=globs)
    print("{}: {} failures out of {} tests".format(modname, fail, total))
    if fail:
        failed_modules.append(mod)

for mod_to_test in MODULES_TO_TEST:
    test_mod(mod_to_test)

sys.exit(1 if failed_modules else 0)
