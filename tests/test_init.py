"""Tests for barcodeforge/__init__.py"""

import unittest
from barcodeforge import __version__


class TestInit(unittest.TestCase):
    def test_example(self):
        self.assertTrue(True)

    def test_version_is_string(self):
        self.assertIsInstance(__version__, str)


if __name__ == "__main__":
    unittest.main()
