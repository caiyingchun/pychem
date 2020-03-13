# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:32:02 2020

@author: 16886
"""


#!/usr/bin/python
#Import required modules
from PySide import QtCore, QtGui
#Import model
from model import SDmodel
#The Molblock viewer class
class MolBlockView(QtGui.QTextBrowser):
    def __init__(self, model, parent=None):
        #Also init the super class
        super(MolBlockView, self).__init__(parent)
        #This sets the window to delete itself when its closed, so it doesn't keep querying the model
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        #We set the model we got on initiation
        self.model = model
        #Set font
        self.font = QtGui.QFont("Monospace")
        self.font.setStyleHint(QtGui.QFont.Monospace)
        self.setFont(self.font)
        self.setGeometry(600, 100, 620, 800)
        #connect the signal to an update of the text
        self.model.selectedChanged.connect(self.updateText)
        #Update text first time
        self.updateText()
    def updateText(self):
        #Function to get the molblock text from the model
        molblock = self.model.getMolBlock()
        self.setText(str(molblock))
if __name__ == "__main__":
    #Import model
    import sys
    from model import SDmodel
    model = SDmodel()
    model.loadSDfile('tester.sdf')
    myApp = QtGui.QApplication(sys.argv)
    molblockview = MolBlockView(model)
    molblockview.show()
    myApp.exec_()