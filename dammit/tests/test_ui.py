#!/usr/bin/env python
from __future__ import print_function

import json
import os
import sys

from unittest import TestCase

from dammit.ui import header, checkbox, paragraph, listing

class TestUIFunctions(TestCase):


    def test_header(self):
        self.assertEquals(header('test'), '# test')
        self.assertEquals(header('test', level=0), '# test')
        self.assertEquals(header('test', level=2), '## test')

    def test_paragraph(self):
        long_text = 'damn' * 40
        self.assertEquals(paragraph(long_text), '\n' + 'damn' * 20 + '\n' +
                                                'damn' * 20 + '\n')

    def test_listing(self):
        d = {'a': 2, 'b': 5, 'c': 10}
        exp_d = '* a: 2\n* b: 5\n* c: 10\n'
        exp_l = '* a\n* b\n* c\n'

        self.assertEquals(listing(d), exp_d)
        self.assertEquals(listing(d.keys()), exp_l)

    def test_checkbox(self):
        self.assertEquals(checkbox('test'), '[ ] test')
        self.assertEquals(checkbox('test', checked=True), '[x] test')
