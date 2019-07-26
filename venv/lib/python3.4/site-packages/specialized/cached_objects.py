#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Specialized cached objects for efficient computation of memory intentive
data structures
"""

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


class CappedDictWithProducer(dict):
    """A specialized subslass of dict with a size cap and a producer function.
    Dictionary grows until cap is reached. Then values are fed from the
    generator function. Before the cap is reached, if a previous unset key is
    requested, again a callback to generator function is used.

    Args:
        producer: a callable :obj: to generate values
        cap: :obj:`int` to set the upper limit of dictionary size
    """
    def __init__(self, producer, cap, *args, **kwargs):
        if not callable(producer):
            raise TypeError(
                'CappedDictWithProducer() arg 1 must be callable object')
        if not isinstance(cap, int):
            raise TypeError(
                'CappedDictWithProducer() arg 2 must be an integer')
        self.producer = producer
        self.cap = cap
        self.count = 0
        super().__init__(*args, **kwargs)

    def __getitem__(self, key):
        try:
            val = super().__getitem__(key)
        except KeyError:
            val = self.producer(key)
            self.__setitem__(key, val)
        return val

    def __setitem__(self, key, val):
        if self.count < self.cap:
            super().__setitem__(key, val)
            self.count += 1

    def __repr__(self):
        dictrepr = super().__repr__()
        return "{}(producer={}, cap={}, count={}, data={})".format(
            type(self).__name__, self.producer, self.cap, self.count, dictrepr)

    def update(self, *args, **kwargs):
        """Implementation of update function"""
        for key, val in dict(*args, **kwargs).items():
            self[key] = val
