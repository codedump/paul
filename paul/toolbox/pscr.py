#!/usr/bin/python

import os

'''
Useful plotscript helpers.
'''

def here(path,base=__file__):
    '''
    Returns a path relative to this script.
    '''
    return os.path.join (os.path.dirname(base), str(path))
