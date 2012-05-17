#!/usr/bin/python

import logging

log = logging.getLogger ('paul')
log.setLevel (logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel (logging.DEBUG)
log.addHandler (ch)

fmt = logging.Formatter('%(asctime)s %(levelname)s(%(name)s): %(message)s')
ch.setFormatter(fmt)

log.debug ("Starting...")

from paul.browser.browserwindow import BrowserWindow
from PyQt4 import QtGui
import sys

testfile = "/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw"

def run():
    app = QtGui.QApplication (sys.argv)
    start_path = "~"
    if (len(sys.argv) > 1 and len(sys.argv[1]) > 0):
        start_path = sys.argv[1]
    win = BrowserWindow(start_path)
    win.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    run()
