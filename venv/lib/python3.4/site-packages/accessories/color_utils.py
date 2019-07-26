# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 11:12:04 2014

@author: barpat
"""

from colour import Color
#import brewer2mpl

try:
  basestring
except NameError:
  basestring = str

PALETTE_TYPES = ('quantitative', 'qualitative', 'sequential', 'diverging', 'unknown')

class Palette(object):
    '''
    A palette class to contain color objects from colour package
    '''
    def __init__(self, name, family, colors=None, pal_type=None):
        if isinstance(name, basestring):
            self.name = name
        else:
            raise TypeError("'name' of the palette has to be a string type")
        if isinstance(family, basestring):
            self.family = family
        else:
            raise TypeError("'family' of palette has to be a string type")
        self.colors = []
        if isinstance(colors, list) or isinstance(colors, tuple):
            self.colors = [colorize(color) for color in colors]
        if pal_type in PALETTE_TYPES:
            self.pal_type = pal_type
        else:
            self.pal_type = 'unknown'

    def __getattr__(self, label):
        if ('get_' + label) in self.__class__.__dict__:
            return getattr(self, 'get_' + label)()
        else:
            raise AttributeError("'%s' not found" % label)

    # getters
    def get_ncols(self):
        return len(self.colors)
    def get_hsl(self):
        return [color.hsl for color in self.colors]
    def get_hex(self):
        return [color.hex for color in self.colors]
    def get_hex_l(self):
        return [color.hex_l for color in self.colors]
    def get_rgb(self):
        return [color.rgb for color in self.colors]

    # pretty
    def __str__(self):
        return "{}".format(self.name)
    def __repr__(self):
        return "<Palette {} ('{}', {} colors)>".format(self.name, self.pal_type, self.ncols)

class color_gradient():
    ''' A class to conveniently return color gradients '''
    valid_gradients = ('linear1', 'linear2', 'bezier')
    valid_returns = ('hex')
    def __init__(self, color_s, gradient_type='linear1', return_type='hex'):
        self.gradient_type = self.check_gradient_type(gradient_type)
        self.return_type = self.check_return_type(return_type)
        if not color_s:
            color_s = ['White', 'Black']
        if isinstance(color_s, tuple):
            color_s = list(color_s)
        if not isinstance(color_s, list):
            color_s = [color_s]
        color_s = [colorize(c) for c in color_s]
        self.edge_colors = color_s
        self.init_color = color_s[0]
        try:
            self.last_color = color_s[-1]
        except IndexError:
            self.last_color = Color('black')
    def check_gradient_type(self, gradient_type):
        if gradient_type and gradient_type in self.valid_gradients:
            return gradient_type
        elif hasattr(self, 'gradient_type'):
            return self.gradient_type
        else:
            return self.valid_gradients[0]
    def check_return_type(self, return_type):
        if return_type and return_type in self.valid_returns:
            return return_type
        elif hasattr(self, 'return_type'):
            return self.return_type
        else:
            return self.valid_returns[0]
    def gradient(self, n=5, gradient_type=None, return_type=None):
        gradient_type = self.check_gradient_type(gradient_type)
        return_type = self.check_return_type(return_type)
        if self.gradient_type == 'linear1':
            cur_gradient = linear_gradient1(self.init_color, self.last_color, n)
        elif self.gradient_type == 'linear2':
            cur_gradient = linear_gradient2(self.init_color, self.last_color, n)
        elif self.gradient_type == 'bezier':
            cur_gradient = bezier_gradient(self.edge_colors, n)
        if self.return_type == 'hex':
            return [c.hex_l for c in cur_gradient]
def colorize(color):
    if isinstance(color, Color):
        return color
    if isinstance(color, int):
        return Color(red=color)
    if isinstance(color, basestring):
        return Color(color)
    if isinstance(color, tuple) or isinstance(color, list):
        if len(color) >= 3:
            return Color(rgb=color[:3])
    raise TypeError("'{}' type can't be colorized.".format(type(color).__name__))

def linear_gradient1(init_color, finish_color=Color('black'), n=10):
    ''' returns a gradient list of (n) colors between two Color objects.
        Based on Ben Southgate's function on http://bsou.io/p/3 '''
    # Starting and ending colors in RGB form
    s = init_color.rgb
    f = finish_color.rgb
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [ s[j] + (float(t)/(n-1))*(f[j]-s[j])
                        for j in range(3)]
        # Add it to our list of output colors
        RGB_list.append(curr_vector)
    return [Color(rgb=rgb) for rgb in RGB_list]

def linear_gradient2(init_color, finish_color=Color('black'), n=10):
    ''' returns a gradient list of (n) colors between two Color objects.
        Based on 'colour' module's internal range_to function '''
    return init_color.range_to(finish_color, n)

fact_cache = {}
def fact(n):
    ''' Memoized factorial function '''
    try:
        return fact_cache[n]
    except(KeyError):
        if n == 1 or n == 0:
            result = 1
        else:
            result = n*fact(n-1)
        fact_cache[n] = result
        return result


def bernstein(t,n,i):
    ''' Bernstein coefficient '''
    binom = fact(n)/float(fact(i)*fact(n - i))
    return binom*((1-t)**(n-i))*(t**i)


def bezier_gradient(colors, num_out=100):
    ''' Returns a "bezier gradient" dictionary using a given list of
        colors as control points. Based on Ben Southgate's function
        on http://bsou.io/p/3'''
    # RGB vectors for each color, use as control points
    RGB_list = [color.rgb for color in colors]
    n = len(RGB_list) - 1

    def bezier_interp(t):
        ''' Define an interpolation function for this specific curve'''
        # List of all summands
        summands = [ list(map(lambda x: bernstein(t,n,i)*x, c))
                     for i, c in enumerate(RGB_list)]
        # Output color
        out = [0.0,0.0,0.0]
        # Add components of each summand together
        for vector in summands:
            for c in range(3):
                out[c] += vector[c]
        return out
    gradient = [bezier_interp(float(t)/(num_out-1)) for t in range(num_out)]
    return [Color(rgb=rgb) for rgb in gradient]
