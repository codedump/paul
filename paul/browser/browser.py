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

def create(path='~', args=[]):
    app = gui.get_app_qt4 (args)

    # find the first non-option argument
    if (len(args) > 1):
        for a in args[1:]:
            if a[0] != '-':
                path = a
                break

    log.debug ("Start path is '%s'" % path)

    main_win = BrowserWindow(path)
    main_win.show()
    if not gui.is_event_loop_running_qt4():
        log.debug ("Starting main event loop")
        app.exec_()
    else:
        log.debug ("Event loop is already running")
    return main_win

if __name__ == "__main__":
    create(args=sys.argv)
