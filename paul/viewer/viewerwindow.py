
import logging
log = logging.getLogger (__name__)

import sys, os, random, imp, subprocess, tempfile
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from paul.loader import igor
from paul.base.wave import Wave
from paul.viewer.matplotlibwidget import MatplotlibWidget
#from paul.viewer.plotscript import *

class PlotscriptToolbar(QtGui.QToolBar):
    '''
    Small widget responsible for selection of a plotscript.
    It contains a dropdown list with several options suitable
    for selecting a plotscript, aswell as a number of
    plotscript-related buttons (Edit, Reload, ...).
    '''

    fileSelected = QtCore.pyqtSignal('QString')

    def __init__ (self, parent=None, name=None, default='(none)'):
        QtGui.QToolBar.__init__ (self, parent)

        self.path_ref  = os.getcwd()
        self.tmp_files = []            # list of temporary files we create --
                                       # should be deleted, but that probably won't work.
        self.file_cur  = ''            # holds the currently selected file name
        self.file_ext  = ".pp"         # plotscript extension
        self.file_hint = os.getcwd()   # holds the reference path (i.e. the path
                                       # that we will use as a hint for different
                                       # locators). can be a directory or a file name,
                                       # is usually the path of the currently displayed
                                       # wave, and is set by the parent widget.

        self.combo = QtGui.QComboBox()
        self.combo.setSizePolicy (QtGui.QSizePolicy.Expanding,
                                  QtGui.QSizePolicy.Minimum)
        self.btn_edit = QtGui.QPushButton("&Edit")
        self.btn_edit.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                                  QtGui.QSizePolicy.Minimum)
        self.btn_edit.clicked.connect (self.onEdit)
        for w in [ self.combo, self.btn_edit ]:
            self.addWidget (w)

        self.combo.addItem ('(none)',         None)
        self.combo.addItem ('Other...',       self.locBrowse)
        self.id_other = 1  # index of the "Other..." entry --
                                               # this one gets special treament.
        self.combo.insertSeparator(2)
        self.combo.addItem ('Automagic',      self.locAuto)
        self.combo.addItem ('File Default',   self.locByFile)
        self.combo.addItem ('Folder Default', self.locByFolder)
        self.combo.addItem ('User Default',   self.locByUser)
        self.combo.insertSeparator(7)
        self.id_user = 8  # index at which we can insert user-selected
                          # plotscripts into the combo box :-)

        self.combo.setCurrentIndex (-1)
        self.combo.currentIndexChanged.connect (self.onSelected)
        def_id = self.combo.findText (default)
        self.combo.setCurrentIndex (def_id) # trigger the default loader

    def __del__(self):
        for i in self.tmp_files:
            log.debug ("Removing temporary plotscript %s" % i)
            os.remove(i)

    #
    # Ps management stuff (resolving script paths, loading modules...)
    #
    def locAuto(self, ref_path):
        '''
        Tries to automagically guess a correct plot script for the current
        wave by trying subsequently:
            . the file plotscript (.pp)
            . the directory plotscript (default.pp in the wave's dir)
            . the directory plotscript of any parent directories
            . the user plotscript
        It returns an empty string if none of the tried paths exist.
        '''
        path = self.locByFile(ref_path)
        if os.path.isfile(path):
            log.debug ("Automagic Ps (file) for %s: %s" % (ref_path, path))
            return path

        base = ref_path
        while True:
            path = self.locByFolder(base)
            if os.path.isfile(path):
                log.debug ("Automagic Ps (folder family) for %s: %s"
                           % (ref_path, path))
                return path
            new_base = os.path.dirname(str(base))
            if new_base == base:
                break
            base = new_base

        path = self.locByUser(ref_path)
        if os.path.isfile(path):
            log.debug ("Automagic Ps (user template) for %s: %s"
                       % (ref_path, path))
            return path

        return ''


    def locByFile(self, ref_path):
        '''
        Return an indivitual plotscript for the currently
        selected wave. The naming convention is <ref_path>.pp.
        If it doesn't exist, return an empty string.
        '''
        path = ref_path+self.file_ext
        log.debug ("File Ps for wave %s: %s" % (ref_path, path))
        return path


    def locByFolder(self, ref_path):
        '''
        Return the default plotscript for the folder in which
        the specified wave is saved. The name of
        the plotscript is 'default.pp'.
        '''
        path = os.path.join (os.path.dirname (str(ref_path)), "default"+self.file_ext)
        log.debug ("Folder Ps for %s: %s" % (ref_path, path))
        return path


    def locByUser(self, ref_path):
        '''
        Return the path of the template plotscript in the current user's
        home directory, or empty string if the script is not available.
        '''
        path = os.path.expanduser ("~/.paul/default"+self.file_ext)
        log.debug ("User plotscript for %s: %s" % (ref_path, path))
        return path


    def locByCombo (self, ref_path):
        '''
        Return the path of the plotscript currently selected by the user.
        The ref_path parameter is ignored. This is the default locator
        script attributed to all user-selected files that appear in the
        combo box. This is why the full path is necessary in the combo box.
        (If we ever decide to show short paths in the combo, we'll need to
        make some kind of mapping for this...).
        '''
        path = str(self.combo.currentText())
        log.debug ("User-selected plotscripts is '%s'" % path)
        return path

    def locBrowse (self, path_hint='default'):
        '''
        Pops up a 'Browse file' dialog box and returns the selected file name.
        '''
        script_file = str(QtGui.QFileDialog.getOpenFileName (self, "Select Paul plot csript", path_hint,
                                                             "Python scripts (*%s)" % self.pscr.file_ext))

        if len(script_file) == 0:  # aborted by user
            return

        if len(script_file) < len(self.file_ext) or script_file[-3:] != self.file_ext:
            script_file = script_file+self.file_ext
        
        if self.combo.findText(script_file) < 0:
            self.combo.insertItem (self.id_user, script_file, self.locByCombo)
                
        return script_file

    
    @QtCore.pyqtSlot()
    def onEdit (self):
        '''
        Called when user clicks the "Edit" button in the viewer window.
        It is supposed to start an editor on the current plotscript file.
        '''
        editor = "emacs"
        if len(self.file_cur) == 0:
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pp", prefix="paul-")
            os.write(tmp_fd,
                     "#!/usr/bin/python\n"
                     "\n"
                     "def decorate (fig, waves):\n"
                     "    pass\n")
            os.close(tmp_fd)
            log.debug ("Loading temporary plotscript %s" % tmp_path)
            self.tmp_files.append (tmp_path)
            self.combo.insertItem (self.id_user, tmp_path, self.locByCombo)
            self.emitPath(tmp_path)
        subprocess.Popen([editor, self.file_cur])
        

    @QtCore.pyqtSlot('int')
    def onSelected (self, key):
        '''
        Called when a plot script is selected from the dropdown box.
        The parameter is the currently selected string. To obtain
        the real path of the script to be loaded, the 'locator'
        has to be called using 'key' as a parameter
        (see psModel() for more information).
        '''

        pscr_path = ''
        locator = self.combo.itemData(key).toPyObject()
        if locator is not None:
            pscr_path = locator(self.path_ref)

        log.debug ("Plot script path is '%s' (locator: %s, combo index: %d)" 
                   % (pscr_path, locator, key))
        self.emitPath(pscr_path)


    def emitPath (self, path):
        '''
        Called when a new file has been selected (or file has been updated).
        This is a wrapper for emitting the corresponding signal and doing
        some related work, to avoid bugs.
        '''
        self.file_cur = path
        self.fileSelected.emit(path)




class ViewerWindow(QtGui.QMainWindow):
    '''
    Class for the main window of the viewer (surprise! :-)
    Displays a plotting canvas and some controls (for plotting,
    plotscript selection etc). Also, plotscripts may attach their
    own toolbars here.
    '''

    # some container classes to get us better organized
    class Plotscript:
        pass
        class Ui:    # container for Plotscript Ui elements
            pass

    class Plotter:
        pass
            

    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100,
                 bgcolor=None, plotscript='Automagic'):
        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle ("Paul Viewer")

        self.pscr = self.Plotscript()
        self.pscr.ui = self.Plotscript.Ui()
        self.plot = self.Plotter()

        self.watcher = QtCore.QFileSystemWatcher()
        self.watcher.fileChanged.connect (self.fileChanged)
        
        self.initMainFrame()
        self.initPlotter()
        self.initScriptLoader(default=plotscript)


    def __del__(self):
        self.cleanup()


    @QtCore.pyqtSlot()
    def cleanup(self):
        '''
        This function is meant to be called when the application / ViewerWindow
        shuts down.
        However, owing to a change of default behavior, there is no way of
        knowing when the window is destroyed: __del__ operator is not called
        anymore when the interpreter shuts down, also prevents the
        destroyed() event from being emitted.
        Until further notice, we'll pack all cleanup tasks in this function,
        but be aware that this still leaves the area un-cleaned...
        '''
        log.debug ("ViewerWindow cleaning up.")
        if hasattr(self.pscr, 'obj') and hasattr(self.pscr.obj, 'exit'):
            self.pscr.obj.exit(self.plot.canvas)


    def initMainFrame(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)
        self.vbox = QtGui.QVBoxLayout()
        self.vbox.setContentsMargins (0, 0, 0, 0)
        self.main_frame.setLayout (self.vbox)


    def initPlotter(self):
        # the plotting area (matplotlib's FigureCanvas and NavigationToolbar)

        self.plot.waves = None        # the currently plotted wave(s)
        self.plot.files = ''
        self.plot.canvas = MatplotlibWidget()
        self.plot.canvas.reset()
        self.plot.tools = NavigationToolbar(self.plot.canvas, self.main_frame)
        self.addToolBar (self.plot.tools)
        for w in [ self.plot.canvas ]:
            self.vbox.addWidget (w)


    def initScriptLoader(self, default='Automagic'):

        self.pscr.cur_file = ''
        self.pscr.toolbar = PlotscriptToolbar(self)
        self.addToolBarBreak()
        self.addToolBar (self.pscr.toolbar)
        self.pscr.toolbar.fileSelected.connect (self.pscrLoad)


    #
    # Matplotlib UI Events
    #
    def mplOnPick(self, event):
        log.debug ("Clicked a bar at %s" 
                   % event.artist.get_bbox().get_points())


    @QtCore.pyqtSlot(Wave)
    def plotWaves (self, wavlist):
        '''
        Plot the specified waves and run the plotscript:
        if 'populate' is available, it will
        take care of plotting the wave(s) itself.
        if only 'decorate' is available, then we are supposed
        to plot the wave(s) and have 'decorate' only do the
        plot decorations afterwards.
        '''
        if wavlist is None or not len(wavlist):
            self.plot.waves = None
            self.setWindowTitle ("Paul Viewer")
            if hasattr(self.plot.canvas, 'axes'):
                self.plot.canvas.reset()
            return

        # If the 'populate' method is available in the plotscript, 
        # we pass all responsibility for plotting to the plotscript.
        if hasattr(self.pscr, 'obj') and hasattr(self.pscr.obj, 'populate'):
                log.debug ("Populating plot (with '%s')." % self.pscr.cur_file)
                self.plot.waves = wavlist
                self.pscr.obj.populate (self.plot.canvas, self.plot.waves)
                self.plot.canvas.draw()

        # ...otherwise we do it ourselves.
        else:
            self.plot.waves = wavlist
            if hasattr(self.pscr, 'obj') and hasattr(self.pscr.obj, 'decorate'):
                self.plot.canvas.reset()
                self.plot.canvas.plot(self.plot.waves, redraw=False)
                log.debug ("Decorating plot (with '%s')." % self.pscr.cur_file)
                self.pscr.obj.decorate (self.plot.canvas, self.plot.waves)
                self.plot.canvas.draw()
            else:
                self.plot.canvas.plot(self.plot.waves, redraw=True)

        if hasattr (self.plot.waves, 'info'):
            self.setWindowTitle ("Paul Viewer: %s" % self.plot.waves.info['name'])
        else:
            self.setWindowTitle ("Paul Viewer <name missing>")


    @QtCore.pyqtSlot('QStringList')
    def plotFiles(self, flist):
        '''
        Loads the specified data file(s) and plots them.
        '''
        if len(flist) < 1:
            log.debug ("Empty list, nothing to do.")
            return

        log.debug ("Plotting %s" % flist[0])
        data = igor.load (flist[0])

        data.info.setdefault('name', os.path.basename(str(flist[0])))
        self.pscr.toolbar.path_ref = str(flist[0])
        self.plot.files = flist
        self.plotWaves (data)

    def refresh(self):
        '''
        This will trigger a plot canvas redraw without all the plotscript
        involved. It is meant to be called when the data changes.
        '''
        self.plot.canvas.draw()


    def replot(self):
        '''
        Wrapper for plot() to easily trigger a replot of the last selection,
        including all the dance (like running the plotscript, etc).
        
        WARNING: This will _not_ reload the files from disk!
                 (Need a reload() method for that?)
        '''
        self.plotWaves(self.plot.waves)


    @QtCore.pyqtSlot('QString')
    def pscrLoad (self, script_file):
        '''
        Loads or re-loads the plotscript, then runs it.

        If the file string is non-zero, and it exists, then it is loaded,
        executed, and registered with the file system watcher.

        If it is non-zero, but non-existing, then just the file system
        watcher is set on it.

        If the file is non-existing, then the current plot script is erased
        from memory.
        
        Called when any of files observed by the QFileSystemWatcher
        is changed on disk. This can be either the plotscript-file,
        or any of the currently displayed wave(s).
        '''

        if len(self.pscr.cur_file) > 0:
            self.watcher.removePath (self.pscr.cur_file)

        if hasattr(self.pscr, 'obj'):
            # tear down the old plotscript
            log.debug ("Killing currently loaded plotscript %s" % script_file)
            self.statusBar().showMessage("Missing plotscript %s" % str(script_file))
            if hasattr(self.pscr, 'obj'):
                if hasattr(self.pscr.obj, 'exit'):
                    log.debug ("exit()'ing old plotscript (%s)" % self.pscr.cur_file)
                    self.pscr.obj.exit(self.plot.canvas, self)
                del self.pscr.obj
                self.pscr.obj = None
            self.plot.canvas.reset()

        if os.path.isfile(str(script_file)):
            # load the script as 'self.plotscript'.
            with open (str(script_file), 'r') as f:                
                log.debug ("Loading '%s'" % str(script_file))
                self.pscr.obj = imp.load_module ('%s.plotscript' % __name__, f, 
                                                  str(script_file), ('', 'r', imp.PY_SOURCE))
                self.statusBar().showMessage ("Plotscript %s" % str(script_file))
                self.watcher.addPath(script_file)
                if hasattr(self.pscr.obj, 'init'):
                    log.debug ("init()'ing plotscript (%s)" % str(script_file))
                    self.pscr.obj.init(self.plot.canvas, self)

        # watch the file (if it's non-zero)
        self.pscr.cur_file = str(script_file)
        self.replot()


    @QtCore.pyqtSlot('QString')
    def fileChanged (self, filename):
        '''
        Called by the file system watcher when the script file has changed.
        '''
        self.pscrLoad (filename)
