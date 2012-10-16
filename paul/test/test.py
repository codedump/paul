#!/usr/bin/python

import paul.loader.igor as igor
from paul.base.wave import *
import numpy as np

testfile1 = "/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw"
testfile2 = "/home/anna/Desktop/ibw/jun09_urs2b_51.ibw"

data, binfo, winfo = igor.loadibw (testfile2)

wav = Wave([1,1])
wav.setScale (0, winfo['sfB'][0], winfo['sfA'][0])
wav.setScale (0, winfo['sfB'][1], winfo['sfA'][1])

print wav

# print data
