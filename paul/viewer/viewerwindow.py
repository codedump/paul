
import logging
log = logging.getLogger (__name__)

import sys, os, random, imp, subprocess, tempfile, glob, shutil
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from paul.loader import igor
from paul.base.wave import Wave
from paul.viewer.matplotlibwidget import MatplotlibWidget


class PlotscriptToolbar(QtGui.QToolBar):
    '''
    Small widget responsible for selection of a plotscript.
    It contains a dropdown list with several options suitable
    for selecting a plotscript, aswell as a number of
    plotscript-related buttons (Edit, Reload, ...).
    '''

    fileSelected = QtCore.pyqtSignal('QString')

    def __init__ (self, parent=None, name=None, default='(none)', pscr_ext=".pp"):
        QtGui.QToolBar.__init__ (self, parent)

        self.path_ref  = os.getcwd()
        self.tmp_files = []            # list of temporary files we create --
                                       # should be deleted, but that probably won't work.
        self.file_cur  = ''            # holds the currently selected file name
        self.file_ext  = pscr_ext      # plotscript extension
        self.file_hint = os.getcwd()   # holds the reference path (i.e. the path
                                       # that we will use as a hint for different
                                       # locators). can be a directory or a file name,
                                       # is usually the path of the currently displayed
                                       # wave, and is set by the parent widget.

        self.combo = QtGui.QComboBox()
        self.combo.setSizePolicy (QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.combo.setMaxVisibleItems (30)
        self.combo.setEditable (True)
        self.combo.setMaximumWidth(300)
        for w in [ self.combo ]:
            self.addWidget (w)
        self.addAction ("Edit", self.onEdit)
        self.addAction ("Reload", self.onForceReload)

        self.combo.addItem ('(none)',         None)
        self.combo.addItem ('Other...',       self.locBrowse)
        self.id_other = self.combo.count()-1   # index of the "Other..." entry --
                                               # this one gets special treament.
        self.combo.insertSeparator(self.combo.count())
        self.combo.addItem ('Automagic',      self.locAuto)
        self.combo.addItem ('File Default',   self.locByFile)
        self.combo.addItem ('Folder Default', self.locByFolder)
        self.combo.addItem ('User Default',   self.locByUser)
        self.combo.insertSeparator(self.combo.count())

        # for convenience, insert the plotscripts from ~/.paul/plugins/
        def_scrpath = os.path.expanduser("~/.paul/plotscripts/*.pp")
        def_scripts = glob.glob(def_scrpath)
        log.debug ("Default plotscripts from %s: %s" % (def_scrpath, def_scripts))
        for f in def_scripts:
            self.combo.addItem (f, self.locByCombo)
        self.combo.insertSeparator(self.combo.count())

        self.id_user = self.combo.count() # index at which we can insert
                                          # user-selected scripts into
                                          # the combo box :-)

        self.combo.setCurrentIndex (-1)
        self.combo.currentIndexChanged.connect (self.onSelected)
        def_id = self.combo.findText (default)
        self.combo.setCurrentIndex (def_id) # trigger the default loader

        # If set to anything else than None, then this widget
        # will be used as a parent for the external editor (if possible).
        self.editor_widget = None
        

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
                                                             "Python scripts (*%s)" % self.file_ext))

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

        ## doesn't work -- bug in emacs ("Unknown command '-parent-id'")
        #if self.editor_widget is None:
        #    # quick hack: create a new frame for emacsclient;
        #    # we want all Paul-related stuff to be opened in the same frame,
        #    # for simplicity...
        #    self.editor_widget = QtGui.QWidget (parent=None)
        
        if self.editor_widget is not None:
            editor = ["emacsclient.emacs-snapshot",  "--parent-id",  "%d" % self.editor_widget.winId()]
        else:
            editor = os.environ.setdefault ('EDITOR', 'gedit').split(" ")
            
            
        if len(self.file_cur) == 0:
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pp", prefix="paul-")
            os.write(tmp_fd,
                     "#!/usr/bin/python\n"
                     "\n"
                     "def decorate (*args, **kwargs):\n"
                     "    '''\n"
                     "    Valid kwargs: can, wav, fig, axes. \n"
                     "    '''\n"
                     "    pass\n")
            os.close(tmp_fd)
            log.debug ("Loading temporary plotscript %s" % tmp_path)
            self.tmp_files.append (tmp_path)
            #self.combo.insertItem (self.id_user, tmp_path, self.locByCombo)
            #self.emitPath (tmp_path)
            self.setPath(tmp_path)
        subprocess.Popen(editor + [self.file_cur])

    @QtCore.pyqtSlot()
    def onForceReload (self):
        '''
        Triggers an unload/reload cycle, by selecting a
        null script (path '') and then the old script
        again. Useful when the old script had syntax errors
        that prevented it from being loaded.
        '''
        old_path = self.file_cur
        self.emitPath('')
        self.emitPath(old_path)


    #@QtCore.pyqtSlot()
    #def onKill (self):
    #    '''
    #    Called when user clicks the "Kill" button in the viewer window.
    #    Selects the empty script, this will tell the parent object
    #    to unconditionally kill the plotscript.
    #    '''
    #    self.emitPath('')


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
        pscr_item = self.combo.itemData(key)
        if pscr_item is None:
            log.error ("Invalid plotscript combo item selected at index %d" % key)
            return

        locator = pscr_item.toPyObject()
        if locator is not None:
            pscr_path = locator(self.path_ref)

        log.debug ("Plot script path is '%s' (locator: %s, combo index: %d)"
                   % (pscr_path, locator, key))
        self.emitPath(pscr_path)


    @QtCore.pyqtSlot('QString')
    def setPath(self, path):
        '''
        Called when a given path is to be visualized as the current
        plot script path. This involves adding the specified path to
        the combo list, and selecting the corresponding item in the list.
        '''
        i = self.combo.findText(path)
        if i < 0:
            self.combo.insertItem (self.id_user, path, self.locByCombo)
            i = self.id_user
        self.combo.setCurrentIndex (i)


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

        class Vars:  # class for 'resident elements', i.e.
            pass     # elements that need to survive between
                     # two reloads of the same plotscript,
                     # but need to be killed on a new plotscript load.

    class Plotter:
        pass


    plotScriptChanged = QtCore.pyqtSignal ('QString')


    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100,
                 bgcolor=None, plotscript='Automagic'):
        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle ("Paul Viewer")

        self.pscr = self.Plotscript()
        self.pscr.ui = self.Plotscript.Ui()
        self.pscr.file_ext = ".pp"
        self.pscr.mod_name = ''
        self.pscr.vars = self.Plotscript.Vars()
        # Check for a variable named 'default_files' in the Plotscript
        # If it exists, and empty list is specified, with plotWaves() or plotFiles(),
        # then load files from that variable for displaying.
        self.pscr.use_default_files = True
        self.plot = self.Plotter()

        # QFileSystemWatcher will send duplicate fileChanged() signals.
        # We do the trick with a QTimer to filter our the bogus signals.
        self.watcher_timer = QtCore.QTimer()
        self.watcher = QtCore.QFileSystemWatcher()
        self.watcher.fileChanged.connect (self.triggerFileChanged)
        self.watcher_timer.setSingleShot(True)
        self.watcher_timer.timeout.connect (self.fileChanged)

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
            self.pscr.obj.exit(can=self.plot.canvas)


    def initMainFrame(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)
        self.vbox = QtGui.QVBoxLayout()
        self.vbox.setContentsMargins (0, 0, 0, 0)
        self.main_frame.setLayout (self.vbox)


    def initPlotter(self):
        # the plotting area (matplotlib's FigureCanvas and NavigationToolbar)

        self.plot.waves = None        # the currently plotted wave(s)
        self.plot.files = []
        self.plot.canvas = MatplotlibWidget()
        self.plot.canvas.reset()
        self.plot.tools = NavigationToolbar(self.plot.canvas, self.main_frame)
        self.plot.tools.setWindowTitle ("Navigation toolbar")
        self.addToolBar (self.plot.tools)
        for w in [ self.plot.canvas ]:
            self.vbox.addWidget (w)


    def initScriptLoader(self, default='Automagic'):

        self.pscr.cur_file = ''
        self.pscr.toolbar = PlotscriptToolbar(self, pscr_ext=self.pscr.file_ext)
        self.pscr.toolbar.setWindowTitle ("Plotscript toolbar")
        self.pscr.toolbar.addAction ("Dump", self.onDump)
        self.addToolBarBreak()
        self.addToolBar (self.pscr.toolbar)
        self.pscr.toolbar.fileSelected.connect (self.loadPlotScript)


    #
    # Matplotlib UI Events
    #
    def mplOnPick(self, event):
        log.debug ("Clicked a bar at %s"
                   % event.artist.get_bbox().get_points())


    @QtCore.pyqtSlot(Wave)
    def plotWaves (self, wavlist, mode='auto'):
        '''
        Plot the specified waves and run the plotscript:
        if 'populate' is available, it will
        take care of plotting the wave(s) itself.
        if only 'decorate' is available, then we are supposed
        to plot the wave(s) and have 'decorate' only do the
        plot decorations afterwards.
        '''
        if wavlist is None or not len(wavlist):
            if hasattr(self.pscr, 'obj') and self.pscr.obj is not None and  hasattr(self.pscr.obj, "default_waves"):
                log.debug ("Plotting <pscr:%s>.default_waves" % os.path.basename(self.pscr.cur_file))
                wavlist = self.pscr.obj.default_waves
            elif hasattr(self.pscr, 'obj') and self.pscr.obj is not None and  hasattr(self.pscr.obj, "default_files"):
                log.debug ("Plotting <pscr:%s>.default_files" % os.path.basename(self.pscr.cur_file))
                flist = [self.pscrRelPath(i) for i in self.pscr.obj.default_files]
                self.plotFiles (flist)
                return
                
        if wavlist is None or not len(wavlist):
            log.info ("Nothing to plot.")
            self.plot.waves = None
            self.setWindowTitle ("Paul Viewer")
            if hasattr(self.plot.canvas, 'axes'):
                self.plot.canvas.reset()
            return

        # If the 'populate' method is available in the plotscript,
        # we pass all responsibility for plotting to the plotscript.
        self.plot.waves = wavlist

        if self.pscrCall ('populate', can=self.plot.canvas,
                          wav=self.plot.waves, fig=self.plot.canvas.fig):
            log.debug ("Plot populated.")
            self.plot.canvas.draw()
        elif self.pscrHasMethod ('decorate'):
            self.plot.canvas.reset()
            self.plot.canvas.plot(self.plot.waves, redraw=False)
            self.pscrCall ('decorate', can=self.plot.canvas, wav=self.plot.waves,
                           fig=self.plot.canvas.fig, axes=self.plot.canvas.axes)
            log.debug ("Plot decorated.")
            self.plot.canvas.draw()
        else:
            log.debug ("Bare plot, no plot script.")
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
        if len(flist) < 1 or (flist is None):
            if hasattr(self.pscr, 'obj') and self.pscr.obj is not None and  hasattr(self.pscr.obj, "default_waves"):
                log.debug ("Plotting <pscr:%s>.default_waves" % os.path.basename(self.pscr.cur_file))
                self.plotWaves (self.pscr.obj.default_waves)
                return
            elif hasattr(self.pscr, 'obj') and self.pscr.obj is not None and  hasattr(self.pscr.obj, "default_files"):
                log.debug ("Plotting <pscr:%s>.default_files" % os.path.basename(self.pscr.cur_file))
                flist = [self.pscrRelPath(i) for i in self.pscr.obj.default_files]
        
                
        if len(flist) < 1 or (flist is None):
            log.debug ("Empty list, nothing to do.")
            return

        log.debug ("File list: %s" % str(flist))
        data = []
        for fname in flist:
            d = igor.load (fname)
            d.info.setdefault('name', os.path.basename(str(fname)))
            self.pscr.toolbar.path_ref = str(fname)
            self.plot.files = list(flist)
            data.append (d)
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


    def pscrRelPath (self, path):
        '''
        Returns a file path relative to the current working directory,
        which was previously assumed to be relative to the plotscript file.
        '''
        if hasattr (self.pscr, 'cur_file'):
            pscr_dir = os.path.dirname(self.pscr.cur_file)
        else:
            pscr_dir = '<unknown>'
            log.warn ('Translation of path (%s) relative to plotscript requested, but no plotscript specified!' % path)
        assert (len(pscr_dir) != 0)
        log.debug ("Path translation for '%s', relative to '%s' requested" % (path, pscr_dir))
        return os.path.normpath(os.path.join(os.path.relpath(pscr_dir, os.getcwd()), path))


    def pscrHasMethod (self, method):
        '''
        Checks if the plotscript has the specified method.
        Returns True or False.
        '''
        return (hasattr(self.pscr, 'obj') and self.pscr.obj is not None and hasattr(self.pscr.obj, method))

        
    def pscrCall (self, method, *args, **kwargs):
        '''
        Runs the specified method in the plotscript. This is just a wrapper
        that implements all due diligence, like checking if the method
        exists (if not, it is ignored) etc.
        Returns 'True' if the method was called (i.e. if it was available),
        or 'False' otherwise.
        Please note that any return value from the method itself is ignored.
        '''
        if self.pscrHasMethod (method):
                proc = getattr(self.pscr.obj, method)
                log.debug ("Calling: '%s' in plotscript '%s'" % (method, self.pscr.cur_file))
                try:
                    proc(*args, **kwargs)
                except:
                    log.error("Failed: '%s' in plotscript '%s'" % (method, self.pscr.cur_file))
                    raise
                return True
        log.debug ("No '%s' in plotscript '%s'" % (method, self.pscr.cur_file))
        return False


    def pscrUnload (self, kill_vars=True):
        '''
        Unloads the current plotscript module. If 'call_exit' is True,
        then the plotscript's exit() method will be called.
        Otherwise it will be ignored.
        '''
        if len(self.pscr.cur_file) > 0:
            self.watcher.removePath (self.pscr.cur_file)

        if hasattr(self.pscr, 'obj') and self.pscr.obj is not None:
            log.debug ("Killing currently loaded plotscript '%s'" % self.pscr.cur_file)
            if kill_vars:
                del self.pscr.vars  # also kill the resident elements
                self.pscr.vars = self.Plotscript.Vars()
            if len(self.pscr.mod_name):
                del sys.modules[self.pscr.mod_name]
            self.pscr.mod_name = ''
            self.pscr.cur_file = ''
            self.pscr.obj = None
            self.plot.canvas.reset()


    @QtCore.pyqtSlot('QString')
    def loadPlotScript(self, script_file):
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

        is_reload = (os.path.abspath(self.pscr.cur_file) == os.path.abspath(str(script_file))) and (self.pscrHasMethod('reload'))

        ## different behavior on reload: don't call init/exit, call only reload
        ## after the new script has been loaded.
        if is_reload:
            self.pscrCall ('unload', can=self.plot.canvas, win=self, fig=self.plot.canvas.fig, vars=self.pscr.vars)
            self.pscrUnload (kill_vars=False)
            self.pscrLoad (str(script_file))
            self.pscrCall ('reload', can=self.plot.canvas, win=self, fig=self.plot.canvas.fig, vars=self.pscr.vars)
        else:
            try:
                self.pscrCall ('exit', can=self.plot.canvas, win=self)
                self.pscrUnload()
                self.pscrLoad(str(script_file))
                self.pscrCall ('init', can=self.plot.canvas, win=self, vars=self.pscr.vars)
                self.plotScriptChanged.emit(script_file)
            except:
                self.pscrUnload()
                raise

        self.replot()  # trigger a replot, this will run the populate/decorate functions


    def setPlotScript (self, path):
        '''
        Sets the specified script as a plotscript.
        This is just a convenience wrapper to be used from the CLI-side
        of things.
        Normally, the plot script is selected by the user, but using this
        function, 'path' is injected at the appropriate position in the
        code path, close to the plotscript-selecting toolbar.
        '''
        log.debug ("Setting plotscript '%s'" % str(path))
        self.pscr.toolbar.setPath (path)


    def pscrLoad (self, sfn):
        '''
        Loads the plotscript as the specified module. Does not execute the init() function.
        The module name is uniquely generated for this instance of the ViewerWindow.
        '''
        if os.path.isfile(sfn):
            with open (sfn, 'r') as f:
                mod_name = '%s.%d_%s' % (__name__, id(self), os.path.basename(sfn[:sfn.rfind(self.pscr.file_ext)]))
                log.debug ("Loading '%s' as module '%s'" % (sfn, mod_name))
                self.pscr.obj = imp.load_module (mod_name, f, sfn, ('', 'r', imp.PY_SOURCE))
                self.pscr.mod_name = mod_name
                self.pscr.cur_file = sfn
                self.watcher.addPath(sfn)

    @QtCore.pyqtSlot('QString')
    def triggerFileChanged(self, filename):
        '''
        Called when a file on the file system has changed.
        Save the file name and trigger a timer signal to
        process it. We do the trick with the timer to avoid
        multiple calls.
        '''
        log.debug ("File changed on disk: %s" % filename)
        if str(filename) == str(self.pscr.cur_file):
            self.watcher_reload_pscr = True
            self.watcher_timer.start(100)
        else:
            log.error ("Unknown file changed: %s" % filename)


    @QtCore.pyqtSlot('QString')
    @QtCore.pyqtSlot()
    def fileChanged (self, filename=''):
        '''
        Called by the file system watcher when the script file has changed.
        '''
        if len(str(filename)) == 0 and self.watcher_reload_pscr == True:
            self.watcher_reload_pscr = False
            filename = self.pscr.cur_file
            log.debug ("Reloading plotscript file: '%s'" % filename)
            self.loadPlotScript (filename)
        else:
            log.error ("Not supposed to be here...")


    @QtCore.pyqtSlot()
    @QtCore.pyqtSlot('QString')
    def onDump(self, path=None):
        '''
        Called when the "Dump" button on the toolbar is pressed. Creates a copy
        of the currently selected plotscript at a user-specified location (pops
        up a file select box to choose the location), and appends the current
        wave files to the file as the 'default_files' variable.

        The purpose is a poor man's 'Save Figure' replacement that will save
        not the image file, but the actual commands to generate it :-)
        '''

        if path is None:
            save_path = str(QtGui.QFileDialog.getSaveFileName (self, "Save plotscript copy as...", os.curdir,
                                                               "Python scripts (*%s)" % self.pscr.file_ext))
        else:
            save_path = str(path)

        if len(save_path) == 0:
            log.debug ("Aborted" % save_path)
            return

        log.debug ("Dumping current plotscript to '%s'" % save_path)

        shutil.copyfile (self.pscr.cur_file, save_path)

        new_paths = [os.path.relpath(str(i), os.path.dirname(save_path)) for i in self.plot.files]
        
        sf = open (save_path, "a")
        sf.write ("default_files = %s\n" % str(new_paths))
        sf.close()

        
    def figSave (file):
        '''
        Poor man's "save figure" implementation: API wrapper around
        the plotscript dump functionality.
        '''
        return self.onDump (file)
