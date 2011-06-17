
import sys, os, random
from qt import *

import matplotlib
matplotlib.use('Agg')
import pylab

from matplotlib.numerix import arange, sin, pi
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from paul.loader import igor


class MatplotlibWidget(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100, bgcolor=None):
	self.parent = parent
	if self.parent:
		bgc = parent.backgroundBrush().color()
		bgcolor = float(bgc.red())/255.0, float(bgc.green())/255.0, float(bgc.blue())/255.0

        self.fig = Figure(figsize=(width, height), dpi=dpi,
                          facecolor=bgcolor, edgecolor=bgcolor)
        self.axes = self.fig.add_subplot(111)

        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()
        self.plot("/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw")
        
        FigureCanvas.__init__(self, self.fig)
        self.reparent(parent, QPoint(0, 0))

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding,
                                         QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

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
                             extent=[dim1_start, dim1_end,
                                     dim2_start, dim2_end])
        else:
            print "MatplotWidget::plot: not implemented for dim > 2"
