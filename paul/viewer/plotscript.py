import logging, sys, os
from PyQt4 import QtCore, QtGui


log = logging.getLogger(__name__)
        

class PlotscriptModel(QtCore.QAbstractItemModel):
    '''
    Obsolete.
    '''
    def __init__ (self, parent = 0):
        QtCore.QAbstractItemModel.__init__(self)
        self.script_list = []
        
        self.addScript ("File Default")
        self.addScript ("User Default")
        self.addScript ("Folder Default")
        self.addScript ("Family Default")

    def addScript (self, key, loc=None):
        '''
        Adds a new script to the script list.
        'locator'  is an object that, if called with the 'descr'-argument,
                   if will generate the full file path of the script file.
                   It needs to understand
        'key'      Is the key (i.e. "name") of the plot script. In most
                   cases, this will be the full path of the script,
                   and the 'locator' will simply return this string.
                   But there are some special keys, like:

                     "File Default"
                         representing the plotscript for the currently
                         selected (displayed) wave, if available

                     "Directory default"
                         the default plotscript for the directory in which
                         the currently selected wave resides
                       
                     "Family default"
                         the 'Directory default' of the currently selected wave,
                         if available, or the 'Directory default' of any parent
                         directory, if not

                     "Template Default"
                         the ~/.paul/plotscript_template.py file.
        '''
        self.script_list.append ({'key': key,
                                  'loc': loc})

    def columnCount (self, parent):
        return 1

    def rowCount (self, parent):
        return len(self.script_list)

    def index(self, row, col, parent):
        # FIXME: Directly indexing a list member... *ugh* this is very ugly :-(
        if not parent.isValid():
            return self.createIndex (row, col, self.script_list[row]['key'])
        else:
            return QtCore.QModelIndex()

    def parent (self, child_index):
        # Only top-level items are present (i.e. no parents), so we return an invalid index.
        return QtCore.QModelIndex()

    def data (self, index, role):
        if (role == QtCore.Qt.DisplayRole):
            return index.internalPointer()
        else:
            return None
