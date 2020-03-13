# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:31:30 2020

@author: 16886
"""


#!/usr/bin/python
#SDbrowser viewer
import sys, time
from PySide.QtGui import *
from PySide.QtCore import QByteArray
from PySide import QtSvg
#Import model
from model import SDmodel
from molblockview import MolBlockView
#The Main browser windows
class MainWindow(QMainWindow):
    def __init__(self,  fileName=None):
        super(MainWindow,self).__init__()
        self.fileName = fileName
        self.filter = "SD files (*.sdf *.sd)"
        self.model = SDmodel()
        #Setup the user interface
        self.initUI()
        #If we get a filename, load it into the model
        if self.fileName != None:
            self.model.loadSDfile(fileName)
    #Update the central widget with SVG from the model
    def update_mol(self):
        self.center.load(QByteArray(self.model.getMolSvg()))
    #Open a new file
    def openFile(self):
        self.fileName, self.filter = QFileDialog.getOpenFileName(self, filter=self.filter)
        self.model.loadSDfile(str(self.fileName))
    #Increment the selected mol with 1
    def nextMol(self):
        #Increment the selected molecule in the model by 1
        self.model.setSelected(self.model.selected + 1)
    #Decrement the selected mol with 1
    def prevMol(self):
        #Decrement the selected molecule in the model by 1
        self.model.setSelected(self.model.selected - 1)
    #Launch the molblockviewer
    def viewMolBlock(self):
        self.molblockBrowser = MolBlockView(self.model)
        self.molblockBrowser.show()
    #Setup the user interface
    def initUI(self):
        #Set Window properties
        self.setWindowTitle("A Simple SD file browser")
        self.setWindowIcon(QIcon('Peptide.png'))
        self.setGeometry(100, 100, 200, 150)
        #Set Central Widget
        self.center = QtSvg.QSvgWidget()
        self.center.setFixedSize(350,350)
        self.setCentralWidget(self.center)
        #Setup the statusbar
        self.myStatusBar = QStatusBar()
        #A permanent widget is right aligned
        self.molcounter = QLabel("-/-")
        self.myStatusBar.addPermanentWidget(self.molcounter, 0)
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage('Ready', 10000)
        #Make the Actions
        self.openAction = QAction( QIcon('Open Folder.png'), 'O&pen',
                                  self, shortcut=QKeySequence.Open,
                                  statusTip="Open an SD file",
                                  triggered=self.openFile)
        self.molblockAction = QAction( QIcon('Page Overview 3.png'), 'V&iew MolBlock',
                                  self, shortcut="Ctrl+M",
                                  statusTip="View MolBlock",
                                  triggered=self.viewMolBlock)
        self.exitAction = QAction( QIcon('Exit.png'), 'E&xit',
                                   self, shortcut="Ctrl+Q",
                                   statusTip="Close the Application",
                                   triggered=self.exit)
        self.prevAction = QAction( QIcon('Left.png'),'Previous', self,
                                   shortcut=QKeySequence.MoveToPreviousChar,
                                   statusTip="Previous molecule",
                                   triggered=self.prevMol)
        self.nextAction = QAction( QIcon('Right.png'),'Next', self,
                                   shortcut=QKeySequence.MoveToNextChar,
                                   statusTip="Next molecule",
                                   triggered=self.nextMol)
        self.aboutAction = QAction( QIcon('Info.png'), 'A&;bout',
                                    self, statusTip="Displays info about SDbrowser",
                                   triggered=self.aboutHelp)
        self.aboutQtAction = QAction("About &Qt", self,
                                statusTip="Qt library About box",
                                triggered=qApp.aboutQt)
        #Setup the menu
        self.fileMenu = self.menuBar().addMenu("&File")
        self.helpMenu = self.menuBar().addMenu("&Help")
        #Setup the Toolbar
        self.mainToolBar = self.addToolBar('Main')
        #Populate the Menu with Actions
        self.fileMenu.addAction(self.openAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)
        self.helpMenu.addAction(self.aboutAction)
        self.helpMenu.addSeparator()
        self.helpMenu.addAction(self.aboutQtAction)
        #Populate the Toolbar with actions.
        self.mainToolBar.addAction(self.openAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.prevAction)
        self.mainToolBar.addAction(self.nextAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.molblockAction)
        #Connect model signals to UI slots
        #Update central widget if the selected molecule changes
        self.model.selectedChanged.connect(self.update_mol)
        #Update the permanent widget in the status bar, if status changes
        self.model.statusChanged.connect(self.molcounter.setText)
        #Finally! Show the UI!
        self.show()
    def exit(self):
        response = QMessageBox.question(self,"Confirmation","This will exit the SD browser\nDo you want to Continue",
                                        QMessageBox.Yes | QMessageBox.No)
        if response == QMessageBox.Yes:
            sdBrowser.quit()
        else:
            pass
    def aboutHelp(self):
        QMessageBox.about(self, "A Basic SD browser",
                "A Simple SD browser where you can see molecules\nIcons from icons8.com")
if __name__ == '__main__':
    # Exception Handling
    try:
        sdBrowser = QApplication(sys.argv)
        #Load with file if provided
        if len(sys.argv) > 1:
            mainWindow = MainWindow(fileName = sys.argv[1])
        else:
            mainWindow = MainWindow()
        sdBrowser.exec_()
        sys.exit(0)
    #Basic Exception handling
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing")
    except Exception:
        print(sys.exc_info()[1])