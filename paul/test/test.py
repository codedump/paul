#!/usr/bin/python

import paul.loader.igor as igor
import numpy as np

testfile = "/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw"

data = igor.loadibw (testfile)

print data
