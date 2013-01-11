#!/usr/bin/python

import logging
import IPython.lib.guisupport as gui
from IPython.frontend.terminal.embed import InteractiveShellEmbed
import IPython.lib.inputhook

from paul.browser.browserwindow import BrowserWindow
from PyQt4 import QtGui
import sys

def create(path='~', args=[], log=None, shell=False):
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

    if log is not None:
        log.debug ("Start path is '%s'" % path)

    main_win = BrowserWindow(path)
    main_win.show()

    if shell:
        # Start an interactive shell, this will be the name space
        # of the shell. Basically, everything needs to be accessed
        # through main_win (aliased 'P' in the shell :-) )

        import numpy as np
        import scipy as sp
        import matplotlib as mpl
        import paul.loader.igor as igor
        import paul.loader.elmitec as elmitec
        import paul.base.wave as wave
        import paul.base.wave as wave
        import paul.toolbox.arpes as arpes
        import paul.toolbox.arplot as arplot
        import paul.toolbox.atrix as atrix
        import paul.toolbox as tools
        import paul.shell as psh
        
        P = main_win
        ipshell = InteractiveShellEmbed.instance()
        IPython.lib.inputhook.enable_gui (gui='qt')
        ipshell()
        
    if not gui.is_event_loop_running_qt4():
        if log is not None:
            log.debug ("Starting main event loop")
        app.exec_()
    else:
        if log is not None:
            log.debug ("Event loop is already running")
    return main_win

if __name__ == "__main__":

    # start an ipython-shell, if '--shell' is specified
    use_shell = False
    if len(sys.argv) > 1:
        i = sys.argv.index('--shell') 
        if i >= 0:
            use_shell = True
            del sys.argv[i]


    # starting logging system
    log = logging.getLogger ('paul')
    log.setLevel (logging.DEBUG)
    if use_shell:
        ch = logging.FileHandler("/tmp/paul-current.log", "w")
    else:
        ch = logging.StreamHandler()
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: %(module)s.%(funcName)s: %(message)s')
    ch.setFormatter(fmt)
    log.debug ("Starting...")

    create(args=sys.argv, shell=use_shell, log=log)
