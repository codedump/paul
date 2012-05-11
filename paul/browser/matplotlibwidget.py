
import sys, os, random
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from paul.loader import igor

import logging
log = logging.getLogger (__name__)


class MatplotlibWidget(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100, bgcolor=None):
        QtGui.QWidget.__init__(self)
	self.parent = parent

        self.fig = Figure(figsize=(width, height), dpi=dpi, facecolor=bgcolor, edgecolor=bgcolor)
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)         # We want the axes cleared every time plot() is called

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.updateGeometry()

        #self.plotFile("/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw")

    def sizeHint(self):
        w = self.fig.get_figwidth()
        h = self.fig.get_figheight()
        return QtCore.QSize(w, h)

    def minimumSizeHint(self):
        return QtCore.QSize(100, 100)

    def compute_initial_figure(self):
       	pass

    # Called to plot file specified by the given full path
    @QtCore.pyqtSlot ('QString')
    def plotFile (self, filename):
        log.debug ("Plotting %s" % filename)
        data, binfo, winfo = igor.loadibw (filename)
        if (data.ndim == 1):
            log.debug ("1D plot")
            self.axes.plot(data)
        elif (data.ndim == 2):
            log.debug ("2D imshow")
            dim1_start = winfo["sfB"][1]
            dim1_end = dim1_start + winfo["sfA"][1]*data.shape[1]
            dim2_start = winfo["sfB"][0]
            dim2_end = dim2_start + winfo["sfA"][0]*data.shape[0]
            # FIXME: loadibw() inverts the data order on one (or more?) of the dimenstions
            #        This means that our scale comes out upside-down, and possibly
            #        also left-right inverted. Need to fix than in igor.loadibw() !
            self.axes.imshow(data, aspect='auto',
                             extent=[min(dim1_start,dim1_end), max(dim1_start,dim1_end),
                                     min(dim2_start,dim2_end), max(dim2_start,dim2_end)])
            self.draw()
            #self.repaint()
            #self.parent.repaint()
        else:
            log.error ("Not implemented for dim > 2")
