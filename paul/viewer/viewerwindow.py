
import logging
log = logging.getLogger (__name__)

import sys, os, random
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from paul.loader import igor
from paul.viewer.matplotlibwidget import MatplotlibWidget

import imp

class ViewerWindow(QtGui.QMainWindow):
    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100,
                 bgcolor=None):
        QtGui.QWidget.__init__(self, parent)

        self.last_path = "~"
        
        self.initMainFrame()
        self.initPlotter()

        self.setWindowTitle ("Paul Viewer")
        self.plotscript_file = ''
        self.watcher = QtCore.QFileSystemWatcher()
        self.watcher.fileChanged.connect (self.fileChanged)
        self.plot_waves = []   # the currently plotted wave(s)


    def initMainFrame(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)
        self.vbox = QtGui.QVBoxLayout()
        self.main_frame.setLayout (self.vbox)


    def initPlotter(self):
        # the plotting area (matplotlib's FigureCanvas and NavigationToolbar)
        self.plot_canvas = MatplotlibWidget()
        self.plot_tools = NavigationToolbar(self.plot_canvas, self.main_frame)
        for w in [ self.plot_tools, self.plot_canvas ]:
            self.vbox.addWidget (w)

        # the script load/save area (a combo box with a (re)load button)
        self.hbox_plotfile = QtGui.QHBoxLayout()
        self.vbox.addLayout (self.hbox_plotfile)
        

        self.plot_scr_list = QtGui.QComboBox()
        self.plot_scr_list.setSizePolicy (QtGui.QSizePolicy.Expanding,
                                            QtGui.QSizePolicy.Minimum)
        self.plot_scr_load = QtGui.QPushButton("&Load")
        self.plot_scr_load.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                          QtGui.QSizePolicy.Minimum)
        self.plot_scr_load.clicked.connect (self.plotScriptSelect)
        for w in [ self.plot_scr_list, self.plot_scr_load ]:
            self.hbox_plotfile.addWidget (w)

    @QtCore.pyqtSlot('QString')
    def plotFile(self, filename):
        '''
        Loads the specified data file and plots it.
        '''
        log.debug ("Plotting %s" % filename)
        data = igor.load (filename)
        self.plot_canvas.plotWave (data)
        self.setWindowTitle ("Paul Viewer: %s" % os.path.basename(str(filename)))
        self.plot_waves = [data]
        self.plotScriptRun ([data])


    def plotScriptRun(self, waves):
        # if a plotscript is defined, use it to modify the plot
        if hasattr(self, 'plotscript') and not self.plotscript == None:
            if hasattr(self.plotscript, 'prepare'):
                waves = self.plotscript.decorate (self.plot_canvas, waves)
            if hasattr(self.plotscript, 'populate'):
                self.plotscript.populate (self.plot_canvas, waves)
            if hasattr(self.plotscript, 'decorate'):
                self.plotscript.decorate (self.plot_canvas, waves)


    def plotScriptLoad (self, script_file):
        '''
        Loads or re-loads the plotscript.
        '''
        log.debug ("BrowserWindow::plotScriptLoad: Loading %s" % str(script_file))
        try:
            f = open (str(script_file), 'r')
            self.plotscript = imp.load_module ('paul.viewer.plotscript', f, str(script_file), ('', 'r', imp.PY_SOURCE))
            self.statusBar().showMessage ("Plotscript %s" % str(script_file))
            if len(self.plotscript_file) > 0:
                self.watcher.removePath (self.plotscript_file)
            self.plotscript_file = str(script_file)
            self.watcher.addPath (self.plotscript_file)
        finally:
            f.close()


    @QtCore.pyqtSlot()
    def plotScriptSelect (self):
        '''
        Called when user selects the "Load" button in the main window.
        
        Ultimately, it loads a plotting script (i.e. a script that
        will beautify the FigureCanvas we are using).

        (Still fighting whether this will select an input from
        the combo box, or pop up a  file-select window.)
        '''
        script_file = QtGui.QFileDialog.getOpenFileName (self, "Select Paul plot csript",
                                                         self.last_path, "Python scripts (*.py)")
        self.plotScriptLoad (script_file)
        

    @QtCore.pyqtSlot('QString')
    def fileChanged(self, path):
        if str(path) == self.plotscript_file:
            self.plotScriptLoad (self.plotscript_file)
            self.plotScriptRun (self.plot_waves)
