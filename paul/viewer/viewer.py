#!/usr/bin/python

import logging
import IPython.lib.guisupport as gui

from paul.viewer.viewerwindow import ViewerWindow
from paul.base.wave import Wave
import paul.loader.igor as igor
from PyQt4 import QtGui
import sys, os
import numpy

def plot (data, win):
    '''
    Plots the specified data. Data can be specified:
      . as a file
      . as a list of files
      . as a Wave()
      . as a list of Wave()s
      . (ndarrays?)
    '''

    # first, make wavs a list (if it's a single file/wave)...
    wavs = []
    if not type(data) is list:
        wavs = [data]
    else:
        wavs = data

    if len(wavs) < 1:
        log.debug ("Empty list, nothing to do.")
        return

    # next, check whether it's a wave or a list of waves.
    if isinstance(wavs[0], numpy.ndarray):
        win.plotWaves (wavs)
    elif isinstance(wavs[0], str):
        ext = wavs[0][wavs[0].rfind("."):]
        
        # check whether we're actually loading a plotscript
        if ext.lower() == ".pp":
            #w = Wave([1, 1])

            # This will trigger reading the variables 'default_waves' or
            # 'default_files' from the Plotscript, if they are available.
            w = None
            win.plotWaves (w)
            win.setPlotScript(os.path.abspath(wavs[0]))
        else:
            # interpret path as a regular IBW path
            win.plotFiles (wavs)
    else:
        log.error ("Don't know what to do with %s." % wavs)


def create (data=[], args=[], log=None):
    '''
    Starts a viewer instance.
    '''
    app = gui.get_app_qt4 (args)
    main_win = ViewerWindow()
    plot (data, win=main_win)
    main_win.show()
        
    if not gui.is_event_loop_running_qt4():
        app.exec_()
        
    return main_win


if __name__ == "__main__":

    log = logging.getLogger ('paul')
    log.setLevel (logging.DEBUG)
    
    ch = logging.StreamHandler()
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)

    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: %(module)s.%(funcName)s: %(message)s')
    ch.setFormatter(fmt)

    log.debug ("Starting...")

    
    if len(sys.argv) <= 1:
        log.error ("Usage: %s <file.ibw>" % sys.argv[0])
        create(args=sys.argv)
    else:
        create(sys.argv[1:], args=sys.argv)
