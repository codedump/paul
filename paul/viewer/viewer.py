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

from paul.viewer.viewerwindow import ViewerWindow
from PyQt4 import QtGui
import sys

canvas = None

def run():
    if len(sys.argv) <= 1:
        log.error ("Usage: %s <file.ibw>" % sys.argv[0])
        sys.exit(-1)

    app = QtGui.QApplication (sys.argv)
    main_win = ViewerWindow()
    main_win.plotFile (sys.argv[1])
    main_win.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    run()
