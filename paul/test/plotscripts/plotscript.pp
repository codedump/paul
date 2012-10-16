#!/usr/bin/python

#
# This is a test plot script. It may contain the following procedures:
#


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
    fig.axes.set_title ('Example Plotscript...')
    fig.draw()
