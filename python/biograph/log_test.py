# pylint: disable=missing-docstring

from __future__ import print_function

import unittest
import logging
from biograph import log_build_stamp


class SaveHandler(logging.Handler):

    def __init__(self):
        logging.Handler.__init__(self)
        self.msgs = []

    def emit(self, msg):
        self.msgs.append(msg.msg)


class LogTestCases(unittest.TestCase):

    def setUp(self):
        self.log = logging.getLogger("biograph")
        self.log.setLevel(logging.INFO)
        self.save_msgs = SaveHandler()
        self.log.addHandler(self.save_msgs)

    def tearDown(self):
        self.log.removeHandler(self.save_msgs)

    def test_logging(self):
        log_build_stamp()

        self.assertEqual(len(self.save_msgs.msgs), 1)
        saved_msg = self.save_msgs.msgs[0]

        self.assertTrue(saved_msg.startswith("Built at") or
                        saved_msg.startswith("Unversioned development build"))

if __name__ == '__main__':
    unittest.main(verbosity=2)
