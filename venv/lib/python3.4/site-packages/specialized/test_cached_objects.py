#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""unittest for CappedDict class"""

import unittest
import cached_objects

class TestCappedDictWithProducer(unittest.TestCase):
    """Provides unit tests for CappedDict objects"""
    def setUp(self):
        def sample_producer(ptype='constant'):
            if ptype == 'constant':
                return lambda x: 1
            elif ptype == 'linear':
                return lambda x: 2*x

        self.test_class = cached_objects.CappedDictWithProducer
        self.constant_prod = sample_producer('constant')
        self.linear_prod = sample_producer('linear')
        self.sample_dict = {'a': 5,
                            2: 'b',
                            (1, 2, 4): [1, 2, 4]}


    def test_init_no_args(self):
        with self.assertRaises(TypeError):
            self.test_class()


    def test_init_from_dict(self):
        with self.assertRaises(TypeError):
            self.test_class(self.sample_dict)

    def test_init_proper_no_data(self):
        try:
            self.test_class(self.constant_prod, 10)
        except Exception as e:
            self.fail('{} can not be initialized properly: {}'.format(
                self.test_class.__name__, e))


if __name__ == '__main__':
    unittest.main()
