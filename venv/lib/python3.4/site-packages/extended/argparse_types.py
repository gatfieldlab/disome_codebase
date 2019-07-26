#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extended argument types for argparse
"""

import argparse
import codecs


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


def int_range(i_min, i_max):
    """Defines an integer taht is bounded inclusively by i_min and i_max

    Args:
        i_min: :obj:`int` - mininum value the integer argument can have
        i_max: :obj:`int` - maximum vlaue the integer argument can have
    """
    if not (isinstance(i_min, int) and isinstance(i_max, int)):
        raise TypeError('int_range arguments have to be integers')
    if not i_min <= i_max:
        raise ValueError('int_range arg 1 <= arg 2')

    def _custom_range(arg_int):
        if arg_int.isdigit and int(arg_int) >= i_min and int(arg_int) <= i_max:
            return int(arg_int)
        else:
            raise argparse.ArgumentTypeError(
                'Expected an integer within [{}-{}]'.format(i_min, i_max))
    return _custom_range


def capped_tuple(mins, maxs):
    """Defines a tuple of integers, each limited to mins and maxs by setup.

    Args:
        mins: :obj:`tuple` or :obj:`list` of integers defining the minimum
            value for each integer variable in the argument
        maxs: :obj:`tuple` or :obj:`list` of integers defining the maximum
            value for each integer variable in the argument
    """
    if not (isinstance(mins, (tuple, list)) and
            isinstance(maxs, (tuple, list)) and
            all(isinstance(val, int) for val in list(mins) + list(maxs))):
        raise TypeError(
            'capped_tuple arguments have to be tuple/list of integers')
    if not len(mins) == len(maxs):
        raise ValueError('capped_tuple arguments have to be of same length')

    def _custom_tuple(arg_str):
        try:
            arg_tuple = tuple([int(i) for i in arg_str.split(',')])
            assert len(arg_tuple) == len(mins) == len(maxs)
        except (ValueError, AssertionError):
            raise argparse.ArgumentTypeError(
                'Expected "{}" comma-delimited "integers"'.format(len(mins)))
        try:
            for i, val in enumerate(arg_tuple):
                assert val >= mins[i] and val <= maxs[i]
        except AssertionError:
            raise argparse.ArgumentTypeError(
                '{} is not within limits set by min={} and max={}'.format(
                    arg_tuple, mins, maxs))
        return arg_tuple
    return _custom_tuple


def int_list(arg_str):
    """Defines a list of integers from a comma-delimited string. Used as an
    argparse type

    Args: None when used as an argparse type
    """
    try:
        return [int(i) for i in arg_str.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError('Expected comma-delimited integers')
    except:
        raise argparse.ArgumentTypeError('Could not parse this argument')


def str_list(arg_str):
    """Defines a list of strings from a comma-delimited string. Used as an
    argparse type

    Args: None when used as an argparse type
    """
    return arg_str.split(',')


def mix_int_list(arg_str):
    """Defines a list of integers from a string of comma-delimited mixed list
    of integers and integer ranges. Example: '2;3:6;8' -> [2,3,4,5,6,8]

    Args: None when used an as argparse type
    """
    if arg_str == '':
        return None
    int_ls = []
    comma_words = arg_str.split(',')
    try:
        for comma_word in comma_words:
            col_words = comma_word.split(':')
            if len(col_words) == 1:
                int_ls.append(int(col_words[0]))
            elif len(col_words) == 2:
                int_ls += list(range(int(col_words[0]), int(col_words[1])+1))
            else:
                raise argparse.ArgumentTypeError(
                    'Expected int ranges in form of A:B')
    except ValueError:
        raise argparse.ArgumentTypeError(
            'Expected a mixed list of integers and integer ranges in form of '
            'A,B:C,D,E:F etc')
    return sorted(set(int_ls))


def unescaped_str(arg_str):
    return codecs.decode(str(arg_str), 'unicode_escape')
