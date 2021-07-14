"""
    setup.py

    The test fixtures (ftest_setup and ftest_teardown) are defined here.
"""
from __future__ import print_function

import os
import re
import sys
import time
import shutil
import datetime

from python.functest.utils.core import TEST_RESULTS
from python.functest.utils.defaults import FTEST_OUTPUT

def process_requirements(required):
    """
        Process test/module requirements: none at the moment
    """
    if not required:
        required = {}

    # if 'upload' in required:
    #     for upload in required['upload']:
    #         print 'Uploading required file: %s -> %s' % (upload[0], upload[1])
    #         upload_file(upload[0], upload[1])

    # if 'reference' in required:
    #     for ref in required['reference']:
    #         if not file_exists('/reference/%s' % ref['name']):
    #             print 'Importing reference: %s' % ref['name']
    #             spiral_cmd(['--no-confirm', '--quiet', 'import_reference', ref['path'], ref['name']])

def ftest_module_setup(required=None):
    """
        Process requirements for a given testing module.
    """
    if not required:
        required = {}

    process_requirements(required)

def ftest_setup(test=None, required=None):
    """ Set up the test environment. This is called before each test is executed. """
    process_requirements(required)
    if "PYTHONPATH" in os.environ:
        os.environ["PYTHONPATH"] += ":%s/modules/bindings/python/" % (os.getcwd())
    else:
        os.environ["PYTHONPATH"] = "%s/modules/bindings/python/" % (os.getcwd())
    #if "PYTHONPATH" in os.environ:
    #    os.environ["PYTHONPATH"] += ":%s/python/biograph" % (os.getcwd())
    #else:
    #    os.environ["PYTHONPATH"] = "%s/python/biograph" % (os.getcwd())
    # FTEST_OUTPUT already appends the date and time. If this is a test, test.TEST_OUTPUT
    # should contain the test ID as well.
    if test:
        # Strip folder/function/module info from test.id()
        pattern = re.compile(r'.*\.(\w+TestCases.*)')
        # Every test gets its own subdir under FTEST_OUTPUT
        test.data_dir = '%s/%s' % (FTEST_OUTPUT, pattern.sub(r'\1', test.id()))
        os.makedirs(test.data_dir)
        # By default, remove data_dir if the test succeeds. This can be overridden per-test.
        test.cleanup = True
        # Create a filehandle to a test-specific log file with line buffering
        test.log_fh = open('%s/test.log' % test.data_dir, 'w', 1)
        # A nice log function with timestamping
        def loggit(msg):
            """ Helper for nicer logging """
            current_time = datetime.datetime.now().time()
            if msg:
                return
            test.log_fh.write('%s %s\n' % (current_time.isoformat(), msg))

        test.log = loggit

        test.coverage_log = open('%s/../coverage.log' % test.data_dir, 'a+')

        # Test timing. This is also tracked by nose --xunit, but is less accessible
        test.start_time = time.time()

def ftest_teardown(test=None):
    """ Clean up after testing. This is called just after each test finishes. """
    if test:
        test.log_fh.close()
        test.coverage_log.close()

        status = None
        # Did the test fail or error out?
        exc = sys.exc_info()
        if exc[0]:
            if 'Exception' in str(exc[0]):
                status = 'ERROR'
            else:
                status = 'FAIL'

        else:
            status = 'PASS'
            if test.cleanup:
                shutil.rmtree(path=test.data_dir, ignore_errors=True)

        test.end_time = time.time()
        test.elapsed_time = test.end_time - test.start_time

        TEST_RESULTS[test.id()] = {'status': status, 'time': test.elapsed_time}


if __name__ == '__main__':
    print("setup.py: required library for functional tests.")
