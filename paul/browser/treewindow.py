from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import os.path
import logging
log = logging.getLogger (__name__)

class TreeWindow (QtGui.QMainWindow):

    # Signal emited when user selected one or more files for plotting.
    wavesSelected = QtCore.pyqtSignal('QStringList')

    def __init__ (self, start_path='~'):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")

        # some variables
        self.current_path = ''
        self.selected_paths = []
        self.set_root_pwd = False # if 'True', change into the Root directory

        # create UI elements
        self.initMainFrame()
        self.initToolbar()
        self.initTree()

        # move to default start directory
        self.setRoot (os.path.expanduser(start_path))


    def initMainFrame(self):
        '''Initializes the inner frame of the tree-view window (main browser feature)'''
        self.vbox = QtGui.QVBoxLayout()
        self.vbox.setContentsMargins (0, 0, 0, 0)
        self.main_frame = QtGui.QWidget()
        self.main_frame.setLayout (self.vbox)
        self.setCentralWidget (self.main_frame)


    def initToolbar(self):
        # create general tree browser layout:
        # one giant VBox, holding the buttons (toolbar?) in the upper part
        # and the file tree in the lower part.
        self.tools_hbox = QtGui.QHBoxLayout()
        self.vbox.addLayout (self.tools_hbox)

        self.tools_label = QtGui.QLabel ("Root: ")

        self.tools_treedrop = QtGui.QComboBox ()
        self.tools_treedrop.setSizePolicy (QtGui.QSizePolicy.Expanding,
                                           QtGui.QSizePolicy.Minimum)
        self.tools_treedrop.currentIndexChanged['QString'].connect (self.setRoot)

        self.tools_btnUp = QtGui.QPushButton ("&Up")
        self.tools_btnUp.show()
        self.tools_btnUp.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.tools_btnUp.clicked.connect (self.treeUp)

        #self.tools_btnSh = QtGui.QPushButton ("&Shell")
        #self.tools_btnSh.show()
        #self.tools_btnSh.setSizePolicy (QtGui.QSizePolicy.Minimum,
        #                                QtGui.QSizePolicy.Minimum)
        #self.tools_btnSh.clicked.connect (self.runShell)

        for w in [ self.tools_label, self.tools_treedrop, self.tools_btnUp ]:
            self.tools_hbox.addWidget (w)


    def initTree(self):
        '''Initializes browser specific stuff (tree-view and buttons)'''

        # model for the file system (dir tree)
        self.filesys = QtGui.QFileSystemModel()
        self.filesys.setRootPath(QtCore.QDir.currentPath())
        self.filesys.setFilter (QtCore.QDir.AllDirs | QtCore.QDir.Dirs | QtCore.QDir.Files | QtCore.QDir.NoDotAndDotDot)
        self.filesys.setNameFilters ("*.ibw")
        self.filesys.setNameFilterDisables (False)

        # the file system browser tree
        self.filetree = QtGui.QTreeView()
        self.vbox.addWidget (self.filetree)
        self.filetree.setSelectionMode (QtGui.QAbstractItemView.ExtendedSelection)
        self.filetree.activated.connect(self.itemActivated)
        self.filetree.setModel (self.filesys)
        self.filetree.selectionModel().selectionChanged.connect(self.itemSelected)

    
    @QtCore.pyqtSlot('QString')
    def setRoot (self, path, silent=0):
        '''
        Sets the specified path as the current root path. Typically
        called when the user double-clicked an entry in the tree view.
        If 'silent' is not 0, then the action will be performed silently
        (i.e. no entry will be changed in the path history dropdown).
        '''

        new_path = os.path.abspath(os.path.expandvars(os.path.expanduser(str(path))))
        if (os.path.abspath(new_path) == self.current_path):
            # this usually happens when we call setCurrentIndex() manually,
            # see below...
            return

        #self.statusBar().showMessage (self.current_path)
        index = self.filesys.index (new_path)
        self.filetree.setRootIndex (index)
        self.filetree.scrollTo (index)
        self.current_path = str(self.filesys.filePath(index))

        print "Setting root index:", self.current_path

        if silent != 0:
            return

        # add path to combo box and make the added entry the current one
        combo_index = self.tools_treedrop.findText (self.current_path)
        if combo_index < 0:
            self.tools_treedrop.insertItem (0, self.current_path)
            combo_index = 0
        self.tools_treedrop.setCurrentIndex (combo_index)
        
        if self.set_root_pwd:
            os.chdir (self.current_path)


    def treeMakeDrop (self):
        '''
        Updates the tree dropdown list.
        The tree-dropdown list is a combo box which contains all
        parent nodes from the currently displayed node (the root
        of the tree-view) to the root of the file system,
        for fast selection.
        Additionally, the dropdown may contain recently selected
        list items.
        '''
        pass
        

        
    @QtCore.pyqtSlot()
    def treeUp(self):
        '''
        Moves the root of the tree view one item up (i.e.
        sets the parent of the current root as the new root).
        '''
        new_path = os.path.dirname(self.current_path)
        log.debug ("Moving to parent: %s" % new_path)
        self.setRoot (new_path)

    @QtCore.pyqtSlot()
    def runShell(self):
        '''
        Runs an introspection IPython shell.
        '''
        from IPython import embed
        embed()


    @QtCore.pyqtSlot('QModelIndex')
    def waveActivated (self, index):
        '''
        Called when an item is activated and that item is a
        regular file (i.e. probably a wave).
        If the item is a wave, it is loaded and plotted.
        '''
        finfo = self.filesys.fileInfo(index)
        fpath = self.filesys.filePath(index)
        log.debug ("Activated file: %s" % fpath)


    @QtCore.pyqtSlot('QModelIndex')
    def itemActivated (self, index):
        '''
        Called when an item in the item tree is activated (double-clicked).
        If the item is a wave, it is loaded and plotted. If it's a directory,
        it activated (i.e. descendend into).
        '''
        finfo = self.filesys.fileInfo(index)
        if finfo.isDir():
            self.setRoot (str(self.filesys.filePath(index)))
        if finfo.isFile() and finfo.isReadable():
            self.waveActivated (index)


    # called when the selection changed in the waveList
    @QtCore.pyqtSlot ('QItemSelection', 'QItemSelection')
    def itemSelected (self, sel, unsel):
        '''
        Called when user selects (marks or unmarks) an item in the wave list.
        The general idea is to load and plot the first selected wave (if
        it's 2D), or all selected waves (if they are 1D).
        
        CAUTION: Basically, the current selection is stored as a dictionary,
                 with the file base name (i.e. the part of the path after the
                 last '/') as the dictionary key. This means that
                 two waves with similar names in different directories will
                 misbehave!
        '''

        # add *sel* to selected_paths
        for i in sel.indexes():
            finfo = self.filesys.fileInfo(i)
            fpath = self.filesys.filePath(i)
            fpath_file = os.path.basename (str(fpath))
            fpath_full = os.path.abspath (str(fpath))
            if fpath not in self.selected_paths:
                if finfo.isFile():
                    log.debug ("Selected file: %s" % str(fpath))
                    self.selected_paths.append (fpath_full)
                elif finfo.isDir():
                    log.debug ("Selected dir: %s" % str(fpath))
                else:
                    log.debug ("...what to do with %s?" % str(fpath))

        # remove *unsel* from selected_paths
        for i in unsel.indexes():
            fpath = self.filesys.filePath(i)
            fpath_file = os.path.basename (str(fpath))
            fpath_full = os.path.abspath (str(fpath))
            if fpath_full in self.selected_paths:
                del self.selected_paths[self.selected_paths.index(fpath_full)]

        self.wavesSelected.emit(self.selected_paths)

        # if the column 0 (the one with the names) is too narrow,
        # expand it to fit the names.
        new_width = self.filetree.sizeHintForColumn (0)
        if (self.filetree.columnWidth(0) < new_width):
            self.filetree.setColumnWidth (0, new_width)
