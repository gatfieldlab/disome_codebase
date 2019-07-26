# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:37:35 2013

@author: barpat

Copyright 2013, Bulak Arpat
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys

class ProgressBar:
    def __init__(self, size, bar=40, barchar='='):
        self.size = size
        self.barlength = bar
        self.barchar = barchar
    def start(self):
        self.update(0)
    def updater(self, c, new_line=False):
        self.current = min(c, self.size)
        p = 1.0 * self.current / self.size
        pg = int(p*self.barlength)
        tail = ">"
        eol = ""
        if pg == self.barlength:
            tail = ""
        if new_line:
            tail = ""
            eol = "\n"
        bar = ("[" + self.barchar*pg + tail +
               " "*(self.barlength-len(tail)-pg) + "]")
        sys.stderr.write("\r{:4.0%}{} {:,.0f} {}".format(p, bar, self.current,
                                                         eol))
        sys.stderr.flush()
    def update(self, c):
        self.updater(c)
    def finish(self, c):
        self.updater(c, new_line=True)

if __name__ == "__main__":
    c = 183932840
    a = 0.0
    update_step = 0
    bar_length = 40
    p = ProgressBar(c,50,'o')
    p.start()
    while True:
        a += 1024
        update_step += 1024
        if a > c:
            a = c
            break
        if update_step > 15000:
            p.update(a)
            update_step = 0
    p.finish(a)
