#!/usr/bin/python

from PyQt4 import QtGui

'''
Widget elements typically used in plugins.
'''

class ValueWidget(QtGui.QWidget):
    '''
    Simple container to hold a label and a SpinBox.
    Used for convenience.
    '''
    def __init__ (self, parent=None, label='value', spin=None, vmin=0, vmax=99, vinc=1, vdef=0, slot=None):
        QtGui.QWidget.__init__ (self, parent)
        self.label = QtGui.QLabel (label, self)
        if spin is None:
            spin = QtGui.QSpinBox(self)
        self.spin = spin
        self.spin.setKeyboardTracking (False)  # ...updates are so terribly slow
        self.spin.setAccelerated (True)
        self.spin.setRange (vmin, vmax)
        self.spin.setSingleStep (vinc)
        self.spin.setValue (vdef)
        if slot is not None:
            self.spin.valueChanged.connect (slot)
        self.box = QtGui.QHBoxLayout(self)
        self.box.setContentsMargins (0, 0, 0, 0)
        self.box.setSpacing (0)
        self.box.addWidget (self.label)
        self.box.addWidget (self.spin)

    
    def value(self):
        return self.spin.value()
