#!/usr/bin/python

import os

'''
Useful plotscript helpers.
'''

def here(base=__file__,path):
    '''
    Returns a path relative to this script.
    '''
    return os.path.join (os.path.dirname(base), str(path))
