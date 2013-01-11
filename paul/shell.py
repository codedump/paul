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
import paul.base.wave as wave
import paul.toolbox.arpes as arpes
import paul.toolbox.arplot as arplot
import paul.toolbox.atrix as atrix
import paul.toolbox as tools

import logging, sys


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

#
# default stand-alone behavior: start a browser with
# an IPython shell in the calling terminal.
#
if __name__ == "__main__":

    # starting logging system
    log = logging.getLogger ('paul')
    log.setLevel (logging.DEBUG)
    ch = logging.FileHandler("/tmp/paul-current.log", "w")
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: %(module)s.%(funcName)s: %(message)s')
    ch.setFormatter(fmt)
    log.debug ("Starting...")
    
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = "~"
        
    P = browser.create(path, shell=True, log=log)
