from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
from paul.browser.treewindow import TreeWindow
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import paul.loader.igor as igor

import os.path

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow):
    def __init__ (self, start_path='~'):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")
        self.resize (800, 500)

        # create UI elements
        self.initMainFrame()
        self.initBrowser()
        self.initViewer()
        self.initRootPwd()

        # start synchronizing the current working directory with self.root
        self.rootpwd = True

        # inter-component connections
        self.tree.wavesSelected.connect(self.viewer.plotFiles)

        # move to default start directory
        self.tree.setRoot (start_path)


    def initMainFrame(self):
        self.hbox = QtGui.QHBoxLayout()
        self.hbox.setSpacing (0)
        self.hbox.setContentsMargins(0, 0, 0, 0)
        self.main_frame = QtGui.QWidget()
        self.main_frame.setLayout (self.hbox)
        self.setCentralWidget (self.main_frame)
        self.splitter = QtGui.QSplitter(self.main_frame)
        self.hbox.addWidget (self.splitter)


    def initBrowser(self):
        # model for the file system (dir tree)
        self.tree = TreeWindow()
        self.tree.setParent (self.splitter)


    def initViewer(self):
        self.viewer = ViewerWindow()
        self.viewer.setParent (self.splitter)


    @QtCore.pyqtSlot('QStringList')
    def plotFiles(wfiles):
        # ...
        pass

    def initRootPwd(self):
        '''
        Starts timer for synchronization between self.root and current working directory
        '''
        self.rootpwd_timer = QtCore.QTimer (self)
        self.rootpwd_timer.timeout.connect (self.syncRootPwd)
        self.rootpwd_timer.stop()
        self.rootpwd_timer.setInterval(500)


    @QtCore.pyqtSlot()
    def syncRootPwd(self):
        '''
        Synchronizes self.root with the current working directory, both ways.
        '''
        cur_path = os.getcwd()
        if os.path.abspath(self.root) != os.path.abspath(cur_path):
            self.root = cur_path
        

    #
    # Some useful browser properties
    #

    @property
    def rootpwd(self):
        '''
        Returns 'True' if the browser is in pwd-sync mode, i.e.
        if the current root directory is synced with the current
        working directory of the application.
        '''
        return self.rootpwd_timer.isActive()


    @rootpwd.setter
    def rootpwd(self, val):
        '''
        Starts/stops pwd sync function.
        '''
        if not hasattr(self, 'rootpwd_timer'):
            return

        if val:
            self.rootpwd_timer.start()
            self.tree.set_root_pwd = True
        else:
            self.rootpwd_timer.stop()
            self.tree.set_root_pwd = False


    @property
    def root(self):
        '''
        The root component of the browser.
        '''
        return self.tree.current_path

    @root.setter
    def root(self, path):
        '''
        Sets the root of the browser.
        '''
        return self.tree.setRoot (path)


    def getSel(self, path='.'):
        '''
        Returns the currently selected wave files as a list of paths,
        relative to *path* (default is the current directory). If 
        *path* is None, then the absolute path is returned.
        '''
        if path is None:
            return [ os.path.normpath(os.path.abspath(i)) for i in self.tree.selected_paths ]
        else:
            return [ os.path.normpath(os.path.relpath(i)) for i in self.tree.selected_paths ]


    @property
    def sel(self):
        '''
        Alias for *BrowserWindow.rsel*.
        '''
        return self.getSel(os.curdir)
        

    @property
    def asel(self, path=None):
        '''
        Wrapper for *BrowserWindow.getSel()*, returns absolute file paths.
        '''
        return self.getSel(None)


    @property
    def rsel(self, path=None):
        '''
        Wrapper for *BrowserWindow.getSel()*, returns relative file paths.
        '''
        return self.getSel(os.curdir)

        
    @property
    def ksel(self):
        '''
        The currently selected waves as a dictionary of absolute file paths,
        with key indices being the basename of the file path.
        (Note that files with the same basename in different directories
        will misbehave. Use the *BrowserWindow._sel* property instead.)
        '''
        return dict({os.path.basename(i): os.path.abspath(i) for i in self.getSel(None)})

    @property
    def kwav(self):
        '''
        A dictionary of Wave() objects, using file base names (i.e. 
        the part of the path after the last '/') as keys. It will misbehave
        if different waves in different directories have the same file name.
        For this purpose, use *BrowserWindow._wav* instead.
        '''
        return dict({i: igor.load(j) for i,j in self.ksel})


    @property
    def wav(self):
        '''
        The currently selected waves as a list of Wave() objects.
        Works correctly for similary named waves in different directories.
        '''
        return [igor.load(i) for i in self.sel]
