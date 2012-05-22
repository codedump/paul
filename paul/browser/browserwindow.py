from PyQt4 import QtGui, QtCore

from paul.browser.matplotlibwidget import MatplotlibWidget
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow):
    def __init__ (self, start_path="/home/florin"):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")

        # some defaults
        self.start_path = start_path
        self.last_path = start_path

        # create UI elements
        self.initMainFrame()
        self.initBrowser()
        self.initPlotter()

        # move to default start directory
        self.dirActivated (self.filesys.index(self.start_path))


    def initMainFrame(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)

        self.hbox = QtGui.QHBoxLayout()
        self.main_frame.setLayout (self.hbox)

        self.splitter = QtGui.QSplitter(self.main_frame)
        self.hbox.addWidget (self.splitter)

    def initBrowser(self):
        # model for the file system (dir tree)
        self.filesys = QtGui.QFileSystemModel()
        self.filesys.setRootPath(QtCore.QDir.currentPath())
        self.filesys.setFilter (QtCore.QDir.AllDirs | QtCore.QDir.Dirs | QtCore.QDir.Files | QtCore.QDir.NoDotAndDotDot)
        self.filesys.setNameFilters ("*.ibw")
        self.filesys.setNameFilterDisables (False)

        # the file system browser tree
        self.file_tree = QtGui.QTreeView(self.splitter)
        self.file_tree.activated.connect(self.dirActivated)
        self.file_tree.setModel (self.filesys)

        # list box for data files (waves)
        self.wave_list = QtGui.QListView(self.splitter)
        self.wave_list.setRootIndex (self.filesys.index(QtCore.QDir.currentPath()))
        self.wave_list.setModel (self.filesys)
        self.wave_list.activated.connect(self.waveActivated)
        self.wave_list.selectionModel().selectionChanged.connect(self.waveSelected)


    def initPlotter(self):
        # main layout manager for the plot area
        self.plot_frame = QtGui.QWidget(self.splitter)
        
        self.vbox_plot = QtGui.QVBoxLayout()
        self.plot_frame.setLayout (self.vbox_plot)

        # the plotting area (matplotlib's FigureCanvas and NavigationToolbar)
        self.plot_canvas = MatplotlibWidget()
        self.plot_tools = NavigationToolbar(self.plot_canvas, self.plot_frame)
        for w in [ self.plot_tools, self.plot_canvas ]:
            self.vbox_plot.addWidget (w)

        # the script load/save area (a combo box with a (re)load button)
        self.hbox_plotfile = QtGui.QHBoxLayout()
        self.vbox_plot.addLayout (self.hbox_plotfile)
        

        self.plot_scr_list = QtGui.QComboBox()
        self.plot_scr_list.setSizePolicy (QtGui.QSizePolicy.Expanding,
                                            QtGui.QSizePolicy.Minimum)
        #self.plot_scr_browse = QtGui.QPushButton ("&Browse")
        #self.plot_scr_browse.setSizePolicy (QtGui.QSizePolicy.Minimum,
        #                                    QtGui.QSizePolicy.Minimum)
        self.plot_scr_load = QtGui.QPushButton("&Load")
        self.plot_scr_load.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                          QtGui.QSizePolicy.Minimum)
        self.plot_scr_load.clicked.connect (self.plotScriptLoad)
        for w in [ self.plot_scr_list, self.plot_scr_load ]:
            self.hbox_plotfile.addWidget (w)


    @QtCore.pyqtSlot()
    def plotScriptLoad (self):
        '''
        Called when user selects the "Load" button in the main window.
        
        Ultimately, it loads a plotting script (i.e. a script that
        will beautify the FigureCanvas we are using).

        (Still fighting whether this will select an input from
        the combo box, or pop up a  file-select window.)
        '''
        script_file = QtGui.QFileDialog.getOpenFileName (self, "Select Paul plot csript",
                                                         self.last_path, "Python scripts (*.py)")
        log.debug ("BrowserWindow::plotScriptLoad: Script file is %s" % str(script_file))



    @QtCore.pyqtSlot('QModelIndex')
    def dirActivated (self, index):
        '''
        Called when user selected an entry from the dir list. The corresponding
        directory will be shown in the wave list.
        '''
        self.wave_list.setRootIndex (index)
        self.file_tree.scrollTo (index)
        self.file_tree.resizeColumnToContents(0)

    @QtCore.pyqtSlot('QModelIndex')
    def waveActivated (self, index):
        '''
        Called when an item in the waveList is activated (double-clicked).
        If the item is a wave, it is loaded and plotted. If it's a directory,
        it activated (i.e. descendend into).
        '''
        finfo = self.filesys.fileInfo(index)
        fpath = self.filesys.filePath(index)
        if finfo.isDir():
            self.wave_list.setRootIndex (index)
            self.file_tree.scrollTo (index)
            self.file_tree.resizeColumnToContents(0)
            self.last_path = str(fpath)
        if finfo.isFile() and finfo.isReadable():
            log.info ("Loading %s" % fpath)
            self.plot_canvas.plotFile (fpath)

    # called when the selection changed in the waveList
    @QtCore.pyqtSlot ('QItemSelection', 'QItemSelection')
    def waveSelected (self, sel, unsel):
        '''
        Called when user selects (marks or unmarks) an item in the wave list.
        The general idea is to load and plot the first selected wave (if
        it's 2D), or all selected waves (if they are 1D).
        '''
        if (sel.isEmpty()):
            return
        ind = sel.indexes()[0]
        log.debug ("Selected %s" % self.filesys.filePath(ind))
        if self.filesys.fileInfo(ind).isFile():
            log.debug ("Passing on %s" % self.filesys.filePath(ind))
            self.waveActivated (ind)
