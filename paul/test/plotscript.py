#!/usr/bin/python

#
# This is a test plot script. It may contain the following procedures:
#
#  __init__()   Will be called upon module loading. Can be used, for
#               example, for showing UI elements.
#
#  prepare()    If present, it will be called to prepare a number of
#               waves for plotting.
#
#  populate()   If present, it will be used to plot a number of waves
#               onto a FigureCanvas.
#
#  decorate()   If present, it will be used to decorate a FigureCanvas.
#

import logging

log = logging.getLogger(__name__)

def __init__():
    log.debug ("plotscript: Initialized.")


def decorate (fig, waves):
    '''
    Decorates the specified FigureCanvas.
    Parameters: 'fig'   Is the FigureCanvas object.
                'waves' Is a list with wave (paul.base.wave) objects
                        containing the data. Usually, the first of these
                        waves will be already plotted.
    '''
    log.debug ("plotscript: Decorating '%s'" % 'fo')
    fig.axes.set_title ('Do you like my wave?')
    fig.draw()
