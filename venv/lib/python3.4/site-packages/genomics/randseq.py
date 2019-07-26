#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generates random sequences confirming certain restrictions"""


import sys
import random
from Bio import SeqUtils
from Bio.SeqUtils import MeltingTemp
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

try:
    import RNA
    RNA_IMPORTED = True
except ImportError:
    RNA_IMPORTED = False


STOP_CODONS = set(['TAA', 'TGA', 'TAG', 'UAA', 'UGA', 'UAG'])
DEF_ALPHABET = IUPAC.unambiguous_rna
MIN_LEN = 6
DEF_GC = 50.0
MIN_GC = 40.0
MAX_GC = 60.0


class RandSeqException(Exception):
    """A simple custom Exception class for randseq.py"""
    def __init__(self, m):
        message = "[randseq] {}".format(m)
        super(RandSeqException, self).__init__(message)


class Restriction(object):
    """A simple base Restriction class
    It also implements a Non-empty restriction"""

    @staticmethod
    def _restrict(seq):
        return len(seq) > 0, len(seq)

    def apply(self, seq):
        """Applies the restriction to the seq obj,
        should return a tuple of (True or False, value)"""
        return self._restrict(seq)


class RestrictGC(Restriction):
    """GC content restriction"""
    def __init__(self, min_gc, max_gc):
        if not (isinstance(min_gc, (int, float)) and
                isinstance(max_gc, (int, float))):
            raise RandSeqException(
                'RestrictGC min_gc and max_gc have to be a float or an int')
        if not (min_gc >= 0 and max_gc >= min_gc and max_gc <= 100.0):
            raise RandSeqException(
                'RestrictGC min_gc >=0, max_gc >= min_gc and max_gc <= 100.0')
        self.min_gc = min_gc
        self.max_gc = max_gc

    def _restrict(self, seq):
        gc_content = SeqUtils.GC(seq)
        return (self.min_gc <= gc_content and gc_content <= self.max_gc,
                gc_content)


class RestrictMeltingTemp(Restriction):
    """Melting temperature restriction"""
    def __init__(self, mt_func, min_mt, max_mt):
        if not (isinstance(mt_func, str) and 
                mt_func in ('Wallace', 'GC', 'NN')):
            raise RandSeqException(
                "RestrictMeltingTemp mt_func is one of: 'Wallace', 'GC', 'NN'")
        if not (isinstance(min_mt, (int, float)) and
                isinstance(max_mt, (int, float))):
            raise RandSeqException(
                'RestrictMeltingTemp min_mt & max_mt are either float or int')
        if not max_mt >= min_mt:
            raise RandSeqException(
                'RestrictGC max_mt >= min_mt')
        self.func = getattr(MeltingTemp, 'Tm_'+mt_func)
        self.min_mt = min_mt
        self.max_mt = max_mt

    def _restrict(self, seq):
        melting_temp = self.func(seq)
        return (self.min_mt <= melting_temp and melting_temp <= self.max_mt,
                melting_temp)


class _RestrictPairing(Restriction):
    """Base class for Folding and Dimerization restrictions"""
    def __init__(self, min_free_energy, max_5p_pair=None):
        clsname = self.__class__.__name__
        if not RNA_IMPORTED:
            raise RandSeqException(
                '{} requires proper import of RNA module'.format(clsname))
        if not isinstance(min_free_energy, (int, float)):
            raise RandSeqException(
                '{} free energy arg has to be int or float'.format(clsname))
        if max_5p_pair is not None:
            if not isinstance(max_5p_pair, int):
                raise RandSeqException(
                    '{} max_5p_pair has to be an int'.format(clsname))
            self.max_5p_pair = '(' * max_5p_pair
        else:
            self.max_5p_pair = None
        self.min_free_energy = min_free_energy

    @staticmethod
    def _folder(seq_str):
        return RNA.fold(seq_str)

    def _restrict(self, seq):
        folding_struc, folding_free_energy = self._folder(str(seq))
        pair_5p_pass = (self.max_5p_pair is None or
                        not folding_struc.startswith(self.max_5p_pair))
        return ((folding_free_energy >= self.min_free_energy and pair_5p_pass),
                (folding_struc, folding_free_energy))


class RestrictFolding(_RestrictPairing):
    """Folding (pinhair) free energy restriction"""


class RestrictDimerization(_RestrictPairing):
    """Self-dimerization free energy restriction"""
    @staticmethod
    def _folder(seq_str):
        return RNA.cofold(seq_str + '&' + seq_str)


class RestrictEndSeq(Restriction):
    """Restriction of user-supplied sequence at the 5' end"""
    def __init__(self, terminus, restricted_seq):
        if terminus not in ('5p', '3p'):
            raise RandSeqException(
                'RestrictEndSeq terminus has to be 5p or 3p')
        if not (isinstance(restricted_seq, (list, tuple)) and
                len(restricted_seq) > 0):
            raise RandSeqException(
                'Restrict5pSeq requires a list/tuple of (string) sequences')
        self.lengths = set([len(seq) for seq in restricted_seq])
        self.restricted_seq = restricted_seq
        self.terminus = terminus

    def _restrict(self, seq):
        if self.terminus == '5p':
            seq_subset = [seq[:l] for l in self.lengths]
        else:
            seq_subset = [seq[-l:] for l in self.lengths]
        found = []
        for subseq in seq_subset:
            if subseq in self.restricted_seq:
                found.append(subseq)
        return len(found) == 0, found


class RestrictStops(Restriction):
    """Restriction of any internal stop codons"""
    def __init__(self, frames):
        if isinstance(frames, int):
            self.frames = [frames-1]
        elif isinstance(frames, (tuple, list)):
            try:
                self.frames = [frame-1 for frame in frames]
            except TypeError:
                raise RandSeqException(
                    'RestrictStops frames should be int or list/tuple of ints')
        valid = [frame in range(3) for frame in self.frames]
        if not all(valid):
            raise RandSeqException(
                'RestrictStops frames should to be 1, 2 or 3')

    def _restrict(self, seq):
        found = []
        for frame in self.frames:
            for start in range(frame, len(seq)-frame, 3):
                if seq[start:start+3] in STOP_CODONS:
                    found.append(frame+1)
                    break
        return len(found) == 0, found


class SeqRestrictions(object):
    """Implements various restrictions for sequences"""
    def __init__(self, *args):
        self.restrictions = []
        for arg in args:
            if isinstance(arg, Restriction):
                self.restrictions.append(arg)
            else:
                raise RandSeqException(
                    'SeqRestrictions only accepts instances of Restrictions')

    def apply(self, seq_iter):
        accepted = []
        for seq in seq_iter:
            res = []
            for restriction in self.restrictions:
                res.append(restriction.apply(seq))
            accept, value = zip(*res)
            if all(accept):
                accepted.append((seq, value))
        return accepted


class WalkerRandom:
    """ Walker's alias method for random objects with different probablities
    From http://code.activestate.com/recipes/576564-walkers-alias-method-for-random-objects-with-diffe/
    """

    def __init__(self, weights, keys=None):
        """builds the Walker tables prob and inx for calls to random().
        The weights (a list or tuple or iterable) can be in any order;
        they need not sum to 1.
        """
        n = self.n = len(weights)
        self.keys = keys
        sumw = sum(weights)
        prob = [w * n / sumw for w in weights]  # av 1
        inx = [-1] * n
        short = [j for j, p in enumerate(prob) if p < 1]
        long = [j for j, p in enumerate(prob) if p > 1]
        while short and long:
            j = short.pop()
            k = long[-1]
            # assert prob[j] <= 1 <= prob[k]
            inx[j] = k
            prob[k] -= (1 - prob[j])  # -= residual weight
            if prob[k] < 1:
                short.append(k)
                long.pop()
        self.prob = prob
        self.inx = inx

    def __str__(self):
        """ e.g. "Walkerrandom prob: 0.4 0.8 1 0.8  inx: 3 3 -1 2" """
        probstr = " ".join(["%.2g" % x for x in self.prob])
        inxstr = " ".join(["%.2g" % x for x in self.inx])
        return "WalkerRandom prob: %s  inx: %s" % (probstr, inxstr)

    def random(self):
        """ each call -> a random int or key with the given probability
            fast: 1 randint(), 1 random.uniform(), table lookup
        """
        u = random.uniform(0, 1)
        j = random.randint(0, self.n - 1)  # or low bits of u
        randint = j if u <= self.prob[j] \
            else self.inx[j]
        return self.keys[randint] if self.keys \
            else randint


class RandSeqGenerator(object):
    def __init__(self, alphabet, mean_gc_content=None):
        if not isinstance(alphabet, (IUPAC.IUPACUnambiguousRNA,
                                     IUPAC.IUPACUnambiguousDNA)):
            raise RandSeqException(
                'RandSeqGenerator alphabet has to be IUPACUnambiguousDNA/RNA')
        if mean_gc_content is None:
            mean_gc_content = DEF_GC
        if not isinstance(mean_gc_content, (int, float)):
            raise RandSeqException(
                'RandSeqGenerator mean_gc_content has to be int or float')
        at_weights = (100 - mean_gc_content) / 2
        gc_weights = mean_gc_content / 2
        letters = alphabet.letters
        weights = []
        for letter in letters:
            if letter in 'GC':
                weights.append(gc_weights)
            else:
                weights.append(at_weights)
        self.randwalker = WalkerRandom(weights, letters)
        self.alphabet = alphabet

    def get_seq(self, length):
        if not (isinstance(length, int) and length >= MIN_LEN):
            raise RandSeqException(
                'get_rand_seq 1 arg has to be an int >= {}'.format(MIN_LEN))
        return Seq(''.join([self.randwalker.random() for i in range(length)]),
                   self.alphabet)


def main():
    """Main body of code"""
    print(get_rand_seq(20))


if __name__ == '__main__':
    sys.exit(main())
