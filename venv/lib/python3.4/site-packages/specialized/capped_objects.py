#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class MutableCappedDict(dict):
    def __init__(self, cap, *args, **kwargs):
        if not isinstance(cap, int):
            raise TypeError('CappedDictWithProducer() arg 1 must be an integer')
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
        # print('update', args, kwargs)
        for k, v in dict(*args, **kwargs).items():
            self[k] = v