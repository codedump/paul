#!/usr/bin/python

'''
This is a convenience module to use from an interactive shell
(IPython, Spyder etc). It imports main Paul modules in a common
name space, and also defines some convenience API functions.
'''

import paul.browser.browser as browser
import paul.viewer.viewer as viewer
from paul.base.wave import *
import numpy as np
import scipy as sp
import matplotlib as mpl
import paul.loader.igor as ig

import paul.toolbox as tools


def view(obj=None):
    '''
    Starts a ViewerWindow on *obj*, where *obj* is either a
    file name (wave or plotscript), or a Wave() object.
    '''
    return viewer.create(obj)


def browse(path='~'):
    '''
    Starts a BrowserWindow on *obj*, where *obj* is expected
    to be a starting path.
    '''
    return browser.create(path)

if __name__ == "__main__":
    pass # no stand-alone stuff yet -- maybe start an ipython shell?
