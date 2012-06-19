#!/usr/bin/python

import logging

log = logging.getLogger ('paul')
log.setLevel (logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel (logging.DEBUG)
log.addHandler (ch)

fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: %(module)s.%(funcName)s: %(message)s')
ch.setFormatter(fmt)

log.debug ("Starting...")

from paul.viewer.viewerwindow import ViewerWindow
from PyQt4 import QtGui
import sys
import numpy

canvas = None
main_win = None


def plot (args, win=None):
    '''
    Plots the specified data. Data can be specified:
      . as a file
      . as a list of files
      . as a Wave()
      . as a list of Wave()s
      . (ndarrays?)
    '''

    global main_win
    
    if win == None:
        win = main_win
    
    # first, make wavs a list (if it's a single file/wave)...
    wavs = []
    if not type(args) is list:
        wavs = [args]
    else:
        wavs = args

    # next, check whether it's a wave or a list of waves.
    if isinstance(wavs[0], numpy.ndarray):
        win.plotWaves (wavs)
    elif isinstance(wavs[0], str):
        win.plotFiles(wavs)
    else:
        log.error ("Don't know what to do with %s." % wavs)


def run (argv=[]):
    '''
    Starts a viewer instance.
    '''
    global main_win
    global canvs

    app = QtGui.QApplication (argv)
    main_win = ViewerWindow()
    if len(argv) > 1:
        plot (argv[1:], win=main_win)
    main_win.show()
    canvas = main_win.plot.canvas

    # need to fix this for the new IPyton 0.11 versions...
    # (event loop is started by IPython there).
    sys.exit(app.exec_())


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        log.error ("Usage: %s <file.ibw>" % sys.argv[0])
        sys.exit(-1)
    run (sys.argv)
