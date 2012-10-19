#!/usr/bin/python

import logging
import IPython.lib.guisupport as gui
from IPython.frontend.terminal.embed import InteractiveShellEmbed
import IPython.lib.inputhook

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

def create(path='~', args=[], shell=False):
    '''
    Creates a paul-browser instance, returns the main window object.
    A "paul-browser" is a GUI window featuring a tree-view and a
    paul viewer window. Clicking on an item in the tree view will
    display it in the viewer window.
    Parameters:
      *path*: Starting point of the tree view, if specified
      *args*: Arguments to pass to the paul-browser. You can
              pass qt command line arguments, or a single
              file name, which, if it is an IBW, will be loaded
              and displayed.
      *shell*: If 'True', the paul-browser will feature an interactive
               IPython shell in the main terminal, which can be used
               for introspection and for control of main features
               of the paul-browser. Enjoy! ;-)
               Default is 'False', but the parameter will be set 
               explicitly 'True' when the browser is called as a
               stand-alone application.
    '''
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

    if shell:
        paul = main_win
        ipshell = InteractiveShellEmbed.instance()
        IPython.lib.inputhook.enable_gui (gui='qt')
        ipshell()
        
    if not gui.is_event_loop_running_qt4():
        log.debug ("Starting main event loop")
        app.exec_()
    else:
        log.debug ("Event loop is already running")
    return main_win

if __name__ == "__main__":
    create(args=sys.argv, shell=True)
