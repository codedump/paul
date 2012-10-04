#!/usr/bin/python

import logging
import IPython.lib.guisupport as gui

log = logging.getLogger ('paul')
log.setLevel (logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel (logging.DEBUG)
log.addHandler (ch)

fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: %(module)s.%(funcName)s: %(message)s')
ch.setFormatter(fmt)

log.debug ("Starting...")

from paul.browser.browserwindow import BrowserWindow
from PyQt4 import QtGui
import sys

canvas = None

def run():
    #app = QtGui.QApplication (sys.argv)
    app = gui.get_app_qt4 (sys.argv)
    start_path = "~"
    if (len(sys.argv) > 1 and len(sys.argv[1]) > 0):
        start_path = sys.argv[1]
    main_win = BrowserWindow(start_path)
    main_win.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    run()
