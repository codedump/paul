#!/usr/bin/python

import os
import logging
log = logging.getLogger(__name__)

'''
Useful plotscript helpers.
'''

def here(path,base=__file__):
    '''
    Returns a path relative to this script.
    '''
    return os.path.join (os.path.dirname(base), str(path))


def mpl_main(argv, populate, log=None):
    '''
    Generic __main__ function for Matplotlib based plotscripts.
    This is supposed to make plotsripts independent of the Paul Viewer.
    In this version, the ploscript can be called with additional argument.
    The argument is a file in which to save the plot.
    
    Arguments:
       *argv*:      The application command line
       *populate*:  A callable that accepts a named parameter 'fig',
                    with a matplotlib.figure.Figure object to draw to.
       *log*:       If specified, the logger will be initialized as
                    a stand-alone logger (i.e. a display handler will
                    be initialized, aswell as a decent formatter).

    Accepted arguments for argv:
       Currently, only one argument (the file name) is being used,
       and it is passed directly to fig.savefig() as a single argument.
       Typically, fig.savefig() will guess the file format from the
       extention.
       In the future, more arguments (like DPI?) may be supported.
    '''
    import sys
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    if log is not None:
        log.setLevel (logging.DEBUG)
        ch = logging.StreamHandler()
        fmt = logging.Formatter('%(levelname)s: %(module)s.%(funcName)s: %(message)s')
        ch.setLevel (logging.DEBUG)
        ch.setFormatter(fmt)
        log.addHandler (ch)

    fig = Figure()
    can = FigureCanvasAgg (fig)
    populate(fig=fig)
    
    if len(sys.argv) > 1:
        fig.savefig (sys.argv[1])
    else:
        log.error ("No output file specified!")
