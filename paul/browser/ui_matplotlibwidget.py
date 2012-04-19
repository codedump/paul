
import sys, os, random
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from matplotlib.numerix import arange, sin, pi
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from paul.loader import igor


class MatplotlibWidget(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100, bgcolor=None):
	self.parent = parent

        self.fig = Figure(figsize=(width, height), dpi=dpi, facecolor=bgcolor, edgecolor=bgcolor)
        self.axes = self.fig.add_subplot(111)

        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()
        self.plot("/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw")
        #self.plot (sys.argv[1])
        
        #FigureCanvas.__init__(self, self.fig)
        #FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        #FigureCanvas.updateGeometry(self)

    def sizeHint(self):
        w = self.fig.get_figwidth()
        h = self.fig.get_figheight()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(100, 100)

    def compute_initial_figure(self):
       	t = arange(0.0, 3.0, 0.01)
	s = sin(2*pi*t)
	self.axes.plot(t, s) 

    def plot (self, filename):
        data, binfo, winfo = igor.loadibw (filename)
        if (data.ndim == 1):
            self.axes.plot(data)
        elif (data.ndim == 2):
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
        else:
            print "MatplotWidget::plot: not implemented for dim > 2"
