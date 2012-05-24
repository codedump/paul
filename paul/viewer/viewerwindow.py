
import logging
log = logging.getLogger (__name__)

import sys, os, random
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from paul.loader import igor
from paul.base.wave import Wave
from paul.viewer.matplotlibwidget import MatplotlibWidget
from paul.viewer.plotscript import *

import imp, subprocess

class ViewerWindow(QtGui.QMainWindow):
    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100,
                 bgcolor=None):

        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle ("Paul Viewer")

        self.last_path = "~"
        self.plot_waves = []   # the currently plotted wave(s)
        self.plot_filename = ''
        self.plotscript_file = ''
        self.watcher = QtCore.QFileSystemWatcher()
        self.watcher.fileChanged.connect (self.fileChanged)
        
        self.initMainFrame()
        self.initPlotter()
        self.initScriptLoader()


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


    def initScriptLoader(self):

        # the script load/save area (a combo box with a (re)load button)
        self.hbox_plotfile = QtGui.QHBoxLayout()
        self.vbox.addLayout (self.hbox_plotfile)
        
        self.plot_scr_combo = QtGui.QComboBox()
        self.plot_scr_combo.setSizePolicy (QtGui.QSizePolicy.Expanding,
                                            QtGui.QSizePolicy.Minimum)
        self.plot_scr_load = QtGui.QPushButton("&Edit")
        self.plot_scr_load.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                          QtGui.QSizePolicy.Minimum)
        self.plot_scr_load.clicked.connect (self.plotScriptEdit)
        for w in [ self.plot_scr_combo, self.plot_scr_load ]:
            self.hbox_plotfile.addWidget (w)

            
        self.plot_scr_combo.addItem ('(none)',         None)
        self.plot_scr_combo.addItem ('Other...', self.plotScriptBrowse)
        self.plot_scr_combo_otherindex = 1  # index of the "Other..." entry --
                                            # this one gets special treament.
        self.plot_scr_combo.insertSeparator(2)
        self.plot_scr_combo.addItem ('Automagic',      self.locateAutomagicPlotscript)
        self.plot_scr_combo.addItem ('File Default',   self.locateFilePlotscript)
        self.plot_scr_combo.addItem ('Folder Default', self.locateFolderPlotscript)
        self.plot_scr_combo.addItem ('User Default',   self.locateUserPlotscript)
        self.plot_scr_combo.insertSeparator(7)
        self.plot_scr_combo_userinsert = 8  # index at which we can insert user-selected
                                            # plotscripts into the combo box :-)

        self.plot_scr_combo.setCurrentIndex (-1)
        self.plot_scr_combo.currentIndexChanged.connect (self.plotScriptSelected)
        self.plot_scr_combo.setCurrentIndex (3) # trigger the 'Automagic' loader


    def locateAutomagicPlotscript(self, wave_path):
        '''
        Tries to automagically guess a correct plot script for the current
        wave by trying subsequently:
            . the file plotscript (.pp)
            . the directory plotscript (_default.pp in the wave's dir)
            . the directory plotscript of any parent directories
            . the user plotscript
        It returns an empty string if none of the tried paths exist.
        '''
        path = self.locateFilePlotscript(wave_path)
        if os.path.isfile(path):
            log.debug ("ViewerWindow::locateAutomagicPlotscript:: Automagic Plotscript (file) for %s: %s" % (wave_path, path))
            return path

        base = wave_path
        while True:
            path = self.locateFolderPlotscript(base)
            if os.path.isfile(path):
                log.debug ("ViewerWindow::locateAutomagicPlotscript:: Automagic Plotscript (folder family) for %s: %s" % (wave_path, path))
                return path
            new_base = os.path.dirname(str(base))
            if new_base == base:
                break
            base = new_base

        path = self.locateUserPlotscript(wave_path)
        if os.path.isfile(path):
            log.debug ("ViewerWindow::locateAutomagicPlotscript:: Automagic Plotscript (user template) for %s: %s" % (wave_path, path))
            return path

        return ''


    def locateFilePlotscript(self, wave_path):
        '''
        Return an indivitual plotscript for the currently
        selected wave. The naming convention is <wave_path>.pp.
        If it doesn't exist, return an empty string.
        '''
        path = wave_path+".pp"
        log.debug ("File Plotscript for wave %s: %s" % (wave_path, path))
        return path


    def locateFolderPlotscript(self, wave_path):
        '''
        Return the default plotscript for the folder in which
        the specified wave is saved. The name of
        the plotscript is '_default.pp'.
        '''
        path = os.path.join (os.path.dirname (str(wave_path)), "_default.pp")
        log.debug ("ViewerWindow::locateFolderPlotscript:: Folder Plotscript for %s: %s" % (wave_path, path))
        return path


    def locateUserPlotscript(self, wave_path):
        '''
        Return the path of the template plotscript in the current user's
        home directory, or empty string if the script is not available.
        '''
        path = os.path.expanduser ("~/.paul/_default.pp")
        log.debug ("ViewerWindow::locateUserPlotscript:: User Plotscript for %s: %s" % (wave_path, path))
        return path


    def locateComboPlotscript (self, wave_path):
        '''
        Return the path of the plotscript currently selected by the user.
        The wave_path parameter is ignored.
        '''
        path = str(self.plot_scr_combo.currentText())
        log.debug ("ViewerWindow::locateComboPlotscript:: Previously user-selected Plotscript is '%s'" % path)
        return path
        

    @QtCore.pyqtSlot('QString')
    def plotFile(self, filename):
        '''
        Loads the specified data file and plots it.
        '''
        log.debug ("Plotting %s" % filename)
        data = igor.load (filename)
        self.plot_canvas.plotWave (data)
        self.setWindowTitle ("Paul Viewer: %s" % os.path.basename(str(filename)))
        self.plot_filename = filename
        self.plot_waves = [data]
        self.plotScriptRun ([data])

    def replot(self):
        '''
        Wrapper for plotFile() to easily trigger a replot of the last selection,
        including all the dance (like running the plotscript, etc).
        A the current, it only re-calls plotFile() using the stored file name.
        Later, if the plotting process on its way from file->plot becomes
        more complex, it might do more.
        '''
        if not hasattr(self, 'plot_filename') or not os.path.isfile(self.plot_filename):
            self.plot_canvas.axes.clear()
            self.plot_canvas.draw()
            self.plotScriptRun (Wave([0]))
        else:
            self.plotFile (self.plot_filename)


    def plotScriptRun(self, waves):
        # if a plotscript is defined, use it to modify the plot
        if hasattr(self, 'plotscript') and not self.plotscript == None:
            if hasattr(self.plotscript, 'populate'):
                self.plotscript.populate (self.plot_canvas, waves)
            elif hasattr(self.plotscript, 'decorate'):
                self.plotscript.decorate (self.plot_canvas, waves)
            else:
                log.warn ("Plotscript is useless: has neither 'populate' nor 'decorate' functions.")
        else:
            log.debug ("No plotscript defined")


    def plotScriptLoad (self, script_file):
        '''
        Loads or re-loads the plotscript.
        '''
        log.debug ("ViewerWindow::plotScriptLoad: Loading %s" % str(script_file))
        if not os.path.isfile(str(script_file)):
            # file doesn't exist, so we'll just tear down the plotscript
            self.statusBar().showMessage("Missing plotscript %s" % str(script_file))
            if hasattr(self, 'plotscript'):
                del self.plotscript
        else:
            # ...otherwise we'll load the script as 'self.plotscript'.
            with open (str(script_file), 'r') as f:
                self.plotscript = imp.load_module ('paul.viewer.plotscript.loaded', f, 
                                                   str(script_file), ('', 'r', imp.PY_SOURCE))
                self.statusBar().showMessage ("Plotscript %s" % str(script_file))
                log.debug ("ViewerWindow::plotScriptLoad: Watching '%s'" % str(script_file))
                self.watcher.addPath(script_file)

        # execute this in any case.
        if len(self.plotscript_file) > 0 and self.plotscript_file != str(script_file):
            self.watcher.removePath (self.plotscript_file)
            log.debug ("ViewerWindow::plotScriptLoad: Un-watching '%s'" % self.plotscript_file)
        self.plotscript_file = str(script_file)
        self.plotScriptRun (self.plot_waves)


    @QtCore.pyqtSlot()
    def plotScriptEdit (self):
        '''
        Called when user clicks the "Edit" button in the viewer window.
        It is supposed to start an editor on the current plotscript file.
        '''
        editor = "emacs"
        if len(self.plotscript_file) == 0:
            log.warn ("Nothing to edit!... (do something useful here, like create a new script?)")
            return
        subprocess.Popen([editor, self.plotscript_file])

    @QtCore.pyqtSlot()
    def plotScriptBrowse (self, start_path=''):
        '''
        Called when user selects the "Load" button in the main window.
        
        Ultimately, it loads a plotting script (i.e. a script that
        will beautify the FigureCanvas we are using).

        (Still fighting whether this will select an input from
        the combo box, or pop up a  file-select window.)
        '''
        if len(start_path) == 0:
            start_path = self.last_path
        script_file = QtGui.QFileDialog.getOpenFileName (self, "Select Paul plot csript",
                                                         start_path, "Python scripts (*.pp)")
        if len(str(script_file)) and self.plot_scr_combo.findText (script_file) < 0:
            self.plot_scr_combo.insertItem (self.plot_scr_combo_userinsert,
                                            script_file, self.locateComboPlotscript)

        return str(script_file)
        

    @QtCore.pyqtSlot('QString')
    def fileChanged(self, path):
        '''
        Called when any of files observed by the QFileSystemWatcher
        is changed on disk. This can be either the plotscript-file,
        or any of the currently displayed wave(s).
        '''
        if str(path) == self.plotscript_file:
            self.plotScriptLoad (self.plotscript_file)
        else:
            pass # it's one of the waves...


    @QtCore.pyqtSlot('int')
    def plotScriptSelected (self, key):
        '''
        Called when a plot script is selected from the dropdown box.
        The parameter is the currently selected string. To obtain
        the real path of the script to be loaded, the 'locator'
        has to be called using 'key' as a parameter
        (see PlotscriptModel() for more information).
        '''

        locator = self.plot_scr_combo.itemData(key).toPyObject()
        if locator is None:
            log.debug ("ViewerWindow::plotScriptSelected: Item %d has no locator." % key)
            path = ''
        else:
            log.debug ("ViewerWindow::plotScriptSelected: Item %d has locator %s." % (key, locator))
            path = locator(self.plot_filename)
        log.debug ("ViewerWindow::plotScriptSelected: Plot script path is '%s' (locator: %s, combo index: %d)" 
                   % (path, locator, key))
        self.plotScriptLoad(path)
        self.replot()
