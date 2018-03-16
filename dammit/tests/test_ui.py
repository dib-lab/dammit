# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import json
import os
import sys

from dammit.ui import header, checkbox, paragraph, listing

class TestUIFunctions():


    def test_header(self):
        assert header('test') == '# test'
        assert header('test', level=0) == '# test'
        assert header('test', level=2) == '## test'

    def test_paragraph(self):
        long_text = 'damn' * 40
        assert paragraph(long_text) == \
               '\n' + 'damn' * 20 + '\n' + 'damn' * 20 + '\n'

    def test_listing(self):
        d = {'a': 2, 'b': 5, 'c': 10}
        exp_d = '* a: 2\n* b: 5\n* c: 10\n'
        exp_l = '* a\n* b\n* c\n'

        assert listing(d) == exp_d
        assert listing(d.keys()) == exp_l

    def test_checkbox(self):
        assert checkbox('test') == '- [ ] test'
        assert checkbox('test', checked=True) == '- [x] test'
