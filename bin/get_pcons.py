#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract PCON score from indexed BGZF compressed files

The *pcons file* had to be compressed with bgzip utility and indexed with
index_bgz.py before this script can work on
"""

import sys
import pickle
from Bio import bgzf


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2019, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


class GetPconsException(Exception):
    """A simple custom Exception class for get_pcons.py
    """
    def __init__(self, m):
        self.message = "[get_pcons] {}".format(m)
        super(GetPconsException, self).__init__(self.message)


class FileException(GetPconsException):
    """A simple custom Exception class for get_pcons.py
    """
    pass


class IndexException(GetPconsException):
    """A simple custom Exception class for get_pcons.py
    """
    def __init__(self, uid, pcons_file):
        self.message = "<{}> could not be found in '{}'\n".format(uid,
                                                                  pcons_file)
        super(IndexException, self).__init__(self.message)


class ExtractionException(GetPconsException):
    """A simple custom Exception class for get_pcons.py
    """
    pass


class PconsFile(object):
    def __init__(self, pcons_file):
        try:
            self.pcons_file_handle = bgzf.BgzfReader(pcons_file, 'r',
                                                     max_cache=5000)
        except IOError as e:
            msg = "Could not read Pcons file!\nI/O error({0}): {1}\n".format(
                e.errno, e.strerror)
            raise FileException(msg)
        self.filename = pcons_file
        try:
            index_file = open(pcons_file+'.idx3', 'rb')
            self.uid_index = pickle.load(index_file)
        except IOError as e:
            msg = "Could not read index file!\nI/O error({0}): {1}\n".format(
                e.errno, e.strerror)
            raise FileException(msg)
        except pickle.UnpicklingError as e:
            msg = "Could not unpickle the index file - possibly wrong format!"
            msg += "\nUnpickling error: {}\n".format(e.message)
            raise FileException(msg)
        except Exception as e:
            msg = "Could not read/unpickle the index file - unknown error! : "
            msg += "{}\n".format(e.__repr__())
            raise FileException(msg)        

    def get_pcons(self, uid):
        extract = ''
        try:
            index_pos = self.uid_index[uid]
        except KeyError:
            return (uid, extract, IndexException(uid, self.filename))
        try:
            self.pcons_file_handle.seek(index_pos)
            for line in self.pcons_file_handle:
                cur_id, _ = line.strip().split('\t', 1)
                if uid in cur_id:
                    extract += line
                else:
                    break
            result = (uid, extract, None)
        except Exception as e:
            msg = "Could not extract '{}' indexed at '{}' "
            msg += "- unknown error : {}\n"
            msg = msg.format(uid, index_pos, e)
            result = (uid, extract, ExtractionException(msg))
        return result

    def __str__(self):
        return "PconsFile: {}".format(self.filename)


def main(argv=None):
    """Main function is run when the script is called from CLI
    """
    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
        sys.stderr.write('usage: get_pcons PCONS-FILE UIDs\n')
        return 1
    try:
        pcons_file = PconsFile(argv[1])
    except FileException as e:
        sys.stderr.write(e.message)
        return 1
    try:
        for uid in argv[2:]:
            _, extract, error = pcons_file.get_pcons(uid)
            if error is not None:
                if isinstance(error, IndexException):
                    sys.stderr.write(error.message)
                    continue
                else:
                    raise error
            else:
                sys.stdout.write(extract)
        return 0
    except ExtractionException as e:
        sys.stderr.write(e.message)
        return 1


if __name__ == "__main__":
    sys.exit(main())
