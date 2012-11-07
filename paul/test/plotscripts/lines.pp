#!/usr/bin/python

import matplotlib as mp
from paul.toolbox.mpltrix import *

'''
Plotscript to show a list of 1D waves as lines. Will use the plotwater()
function from paul.toolbox.mpltrix.
'''

def populate (*args, **kwargs):

    wav = kwargs['wav']
    can = kwargs['can']
    can.reset()

    w = [i / i.max() for i in wav]
    lines = plotwater (can.axes, w, offs=(0, 0))

    kwargs['lines'] = lines
    kwargs['ylim'] = (min([i.min() for i in w]), max([i.max() for i in w]))
    decorate (*args, **kwargs)
    

def decorate (*args, **kwargs):

    ax = kwargs['can'].axes
    ax.set_ylim (kwargs['ylim'])
