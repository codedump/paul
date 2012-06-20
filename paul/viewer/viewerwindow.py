
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

class ViewerWindow(QtGui.QMainWindow):

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
        self.watcher.fileChanged.connect (self.pscrOnFileChanged)
        
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
        for i in self.pscr.tmp_files:
            log.debug ("Removing temporary plotscript %s" % i)
            os.remove(i)


    def initMainFrame(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)
        self.vbox = QtGui.QVBoxLayout()
        self.main_frame.setLayout (self.vbox)


    def initPlotter(self):
        # the plotting area (matplotlib's FigureCanvas and NavigationToolbar)

        self.plot.waves = []        # the currently plotted wave(s)
        self.plot.files = ''

        self.plot.canvas = MatplotlibWidget()
        self.plot.tools = NavigationToolbar(self.plot.canvas, self.main_frame)
        for w in [ self.plot.tools, self.plot.canvas ]:
            self.vbox.addWidget (w)


    def initScriptLoader(self, default='Automagic'):

        # the script load/save area (a combo box with a (re)load button)
        self.hbox_plotfile = QtGui.QHBoxLayout()
        self.vbox.addLayout (self.hbox_plotfile)
        self.vbox.setContentsMargins (0, 0, 0, 0)

        self.pscr.obj  = None
        self.pscr.file = ''
        self.pscr.file_ext = ".pp"
        self.pscr.tmp_files = []
        
        self.pscr.ui.combo = QtGui.QComboBox()
        self.pscr.ui.combo.setSizePolicy (QtGui.QSizePolicy.Expanding,
                                                QtGui.QSizePolicy.Minimum)
        self.pscr.ui.btn_edit = QtGui.QPushButton("&Edit")
        self.pscr.ui.btn_edit.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                                  QtGui.QSizePolicy.Minimum)
        self.pscr.ui.btn_edit.clicked.connect (self.pscrOnEdit)
        for w in [ self.pscr.ui.combo, self.pscr.ui.btn_edit ]:
            self.hbox_plotfile.addWidget (w)

            
        self.pscr.ui.combo.addItem ('(none)',         None)
        self.pscr.ui.combo.addItem ('Other...',       self.pscrLocBrowse)
        self.pscr.ui.id_other = 1  # index of the "Other..." entry --
                                               # this one gets special treament.
        self.pscr.ui.combo.insertSeparator(2)
        self.pscr.ui.combo.addItem ('Automagic',      self.pscrLocAuto)
        self.pscr.ui.combo.addItem ('File Default',   self.pscrLocByFile)
        self.pscr.ui.combo.addItem ('Folder Default', self.pscrLocByFolder)
        self.pscr.ui.combo.addItem ('User Default',   self.pscrLocByUser)
        self.pscr.ui.combo.insertSeparator(7)
        self.pscr.ui.id_user = 8  # index at which we can insert user-selected
                                              # plotscripts into the combo box :-)

        self.pscr.ui.combo.setCurrentIndex (-1)
        self.pscr.ui.combo.currentIndexChanged.connect (self.pscrOnSelected)
        def_id = self.pscr.ui.combo.findText (default)
        self.pscr.ui.combo.setCurrentIndex (def_id) # trigger the 'Automagic' loader


    #
    # Matplotlib UI Events
    #

    def mplOnPick(self, event):
        log.debug ("Clicked a bar at %s" 
                   % event.artist.get_bbox().get_points())


    #
    # Ps management stuff (resolving script paths, loading modules...)
    #

    def pscrLocAuto(self, wave_path):
        '''
        Tries to automagically guess a correct plot script for the current
        wave by trying subsequently:
            . the file plotscript (.pp)
            . the directory plotscript (default.pp in the wave's dir)
            . the directory plotscript of any parent directories
            . the user plotscript
        It returns an empty string if none of the tried paths exist.
        '''
        path = self.pscrLocByFile(wave_path)
        if os.path.isfile(path):
            log.debug ("Automagic Ps (file) for %s: %s" % (wave_path, path))
            return path

        base = wave_path
        while True:
            path = self.pscrLocByFolder(base)
            if os.path.isfile(path):
                log.debug ("Automagic Ps (folder family) for %s: %s"
                           % (wave_path, path))
                return path
            new_base = os.path.dirname(str(base))
            if new_base == base:
                break
            base = new_base

        path = self.pscrLocByUser(wave_path)
        if os.path.isfile(path):
            log.debug ("Automagic Ps (user template) for %s: %s"
                       % (wave_path, path))
            return path

        return ''


    def pscrLocByFile(self, wave_path):
        '''
        Return an indivitual plotscript for the currently
        selected wave. The naming convention is <wave_path>.pp.
        If it doesn't exist, return an empty string.
        '''
        path = wave_path+self.pscr.file_ext
        log.debug ("File Ps for wave %s: %s" % (wave_path, path))
        return path


    def pscrLocByFolder(self, wave_path):
        '''
        Return the default plotscript for the folder in which
        the specified wave is saved. The name of
        the plotscript is 'default.pp'.
        '''
        path = os.path.join (os.path.dirname (str(wave_path)), "default"+self.pscr.file_ext)
        log.debug ("Folder Ps for %s: %s" % (wave_path, path))
        return path


    def pscrLocByUser(self, wave_path):
        '''
        Return the path of the template plotscript in the current user's
        home directory, or empty string if the script is not available.
        '''
        path = os.path.expanduser ("~/.paul/default"+self.pscr.file_ext)
        log.debug ("User plotscript for %s: %s" % (wave_path, path))
        return path


    def pscrLocByCombo (self, wave_path):
        '''
        Return the path of the plotscript currently selected by the user.
        The wave_path parameter is ignored. This is the default locator
        script attributed to all user-selected files that appear in the
        combo box. This is why the full path is necessary in the combo box.
        (If we ever decide to show short paths in the combo, we'll need to
        make some kind of mapping for this...).
        '''
        path = str(self.pscr.ui.combo.currentText())
        log.debug ("User-selected plotscripts is '%s'" % path)
        return path


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
        if not len(wavlist):
            self.plot.waves = []
            self.setWindowTitle ("Paul Viewer")
            if hasattr(self.plot.canvas, 'axes'):
                self.plot.canvas.clear()
            return

        # If the 'populate' method is available in the plotscript, 
        # we pass all responsibility for plotting to the plotscript.
        if hasattr(self.pscr, 'obj') and hasattr(self.pscr.obj, 'populate'):
                log.debug ("Populating plot (with '%s')." % self.pscr.file)
                self.plot.waves = wavlist
                self.pscr.obj.populate (self.plot.canvas, self.plot.waves)

        # ...otherwise we do it ourselves.
        else:
            self.plot.waves = wavlist
            self.plot.canvas.plotWave (self.plot.waves[0], redraw=False)
            if hasattr(self.pscr, 'obj') and hasattr(self.pscr.obj, 'decorate'):
                log.debug ("Decorating plot (with '%s')." % self.pscr.file)
                self.pscr.obj.decorate (self.plot.canvas, self.plot.waves)

        self.plot.canvas.draw()
        if hasattr (self.plot.waves[0], 'info'):
            self.setWindowTitle ("Paul Viewer: %s" % self.plot.waves[0].info['name'])
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
        self.plot.files = flist
        self.plotWaves ([data])

    def refresh(self):
        '''
        This will trigger a plot canvas redraw without all the plotscript
        involved. It is meant to be called when the data changes.
        '''
        self.plot.canvas.draw()


    def replot(self):
        '''
        Wrapper for plotWaves() to easily trigger a replot of the last selection,
        including all the dance (like running the plotscript, etc).
        
        WARNING: This will _not_ reload the files from disk!
                 (Need a reload() method for that?)
        '''
        self.plotWaves(self.plot.waves)


    def pscrLocBrowse (self, path_hint='default'):
        '''
        Pops up a 'Browse file' dialog box and returns the selected file name.
        '''
        script_file = str(QtGui.QFileDialog.getOpenFileName (self, "Select Paul plot csript", path_hint,
                                                             "Python scripts (*%s)" % self.pscr.file_ext))

        if len(script_file) == 0:  # aborted by user
            return

        if len(script_file) < len(self.pscr.file_ext) or script_file[-3:] != self.pscr.file_ext:
            script_file = script_file+self.pscr.file_ext
        
        if self.pscr.ui.combo.findText(script_file) < 0:
            self.pscr.ui.combo.insertItem (self.pscr.ui.id_user, script_file, self.pscrLocByCombo)
                
        return script_file


    def pscrRun(self):
        '''
        If a plotscript has been defined, run it on the current plots.
        This is actually just a wrapper for plotWaves(), since
        the plotting of the waves is dependent on which methods the
        plotscript provides.
        '''
        plotWaves(self.plot.waves)


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

        if len(self.pscr.file) > 0:
            self.watcher.removePath (self.pscr.file)

        if hasattr(self.pscr, 'obj'):
            # tear down the old plotscript
            log.debug ("Killing currently loaded plotscript %s" % script_file)
            self.statusBar().showMessage("Missing plotscript %s" % str(script_file))
            if hasattr(self.pscr, 'obj'):
                if hasattr(self.pscr.obj, 'exit'):
                    log.debug ("exit()'ing old plotscript (%s)" % self.pscr.file)
                    self.pscr.obj.exit(self.plot.canvas)
                del self.pscr.obj
                self.pscr.obj = None
            self.plot.canvas.clear()

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
                    self.pscr.obj.init(self.plot.canvas)

        # watch the file (if it's non-zero)
        self.pscr.file = str(script_file)


    @QtCore.pyqtSlot('QString')
    def pscrOnFileChanged (self, filename):
        '''
        Called by the file system watcher when the script file has changed.
        '''
        if str(filename) == self.pscr.file:
            self.pscrLoad(filename)
            self.replot()
        else:
            log.debug ("Don't know what to do with %s." % str(filename))



    @QtCore.pyqtSlot()
    def pscrOnEdit (self):
        '''
        Called when user clicks the "Edit" button in the viewer window.
        It is supposed to start an editor on the current plotscript file.
        '''
        editor = "emacs"
        if len(self.pscr.file) == 0:
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pp", prefix="paul-")
            os.write(tmp_fd,
                     "#!/usr/bin/python\n"
                     "\n"
                     "def decorate (fig, waves):\n"
                     "    pass\n")
            os.close(tmp_fd)
            log.debug ("Loading temporary plotscript %s" % tmp_path)
            self.pscr.tmp_files.append (tmp_path)
            self.pscr.ui.combo.insertItem (self.pscr.ui.id_user, tmp_path, self.pscrLocByCombo)
            self.pscrLoad (tmp_path)
            self.replot()

        subprocess.Popen([editor, self.pscr.file])
        

    @QtCore.pyqtSlot('int')
    def pscrOnSelected (self, key):
        '''
        Called when a plot script is selected from the dropdown box.
        The parameter is the currently selected string. To obtain
        the real path of the script to be loaded, the 'locator'
        has to be called using 'key' as a parameter
        (see psModel() for more information).
        '''

        pscr_path = ''
        wave_path = ''
        locator = self.pscr.ui.combo.itemData(key).toPyObject()

        if len(self.plot.files) > 0:
            wave_path = self.plot.files[0]
        if locator is not None:
            pscr_path = locator(wave_path)

        log.debug ("Plot script path is '%s' (locator: %s, combo index: %d)" 
                   % (pscr_path, locator, key))

        self.pscrLoad(pscr_path)
        self.replot()
