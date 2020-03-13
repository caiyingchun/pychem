# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:29:36 2020

@author: 16886
"""


#QtCore is the nonGUI stuff.
from PySide.QtCore import QObject, Slot
from PySide.QtCore import Signal
#RDKit stuff
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
#The model holding the SDfile data
class SDmodel(QObject): #Inherit from QObject so that Signals can be emitted
    def __init__(self):
        super(SDmodel, self).__init__() #Init the super QObject class so that it works with QT stuff etc.
        self._selected = 0
        self._status = "Ready"
        self.length = 0
        #Set the counter when the selection is changed (non GUI signal slot example)
        self.selectedChanged.connect(self.setCounter)
    #Define a signal that can be emitted when the selected compound is changed, is called pyqtSignal in PyQt
    selectedChanged = Signal(int, name = 'selectedChanged')
    @property
    def selected(self):
        return self._selected
    #This enables us to do model.selected = 2
    @selected.setter
    def selected(self, value):
        self.setSelected(value)
    #this is more easy to use from with PyQt signals
    def setSelected(self, selected):
        #Prevent setting a selected that doesn't exist
        if selected < 0: selected = 0
        if selected > self.length -1: selected = self.length -1
        #Only set the selected if its changed, we get roundtripping otherwise
        if selected != self._selected:
            self._selected = selected
            print "in model: selected set for ", selected
            #Emit the signal that selected has changed
            self.selectedChanged.emit(self._selected)
    #Decorate the function setCounter
    @Slot()
    def setCounter(self):
        self.setStatus('%s/%s'%(self._selected + 1, self.length))
    #A status signal and property
    statusChanged = Signal(str, name = 'statusChanged')
    @property
    def status(self):
        return self._status
    def setStatus(self, status):
        self._status = status
        self.statusChanged.emit(self._status)
    def loadSDfile(self, filename):
        self.filename = filename
        self.SDMolSupplier = Chem.SDMolSupplier(filename)
        self.length = len(self.SDMolSupplier)
        if self.selected == 0:
            self.selectedChanged.emit(self._selected)
        else:
            self.setSelected(0)
    #Better rendering with SVG
    def getMolSvg(self, kekulize=True, calc2Dcoords=True):
        mol = self.SDMolSupplier[self._selected]
        mc = Chem.Mol(mol.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers() or calc2Dcoords:
            rdDepictor.Compute2DCoords(mc)
        drawer = rdMolDraw2D.MolDraw2DSVG(300,300)
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:','')
        return svg
    def getMolBlock(self):
        print "Getting MolBlock"
        return self.SDMolSupplier.GetItemText(self.selected)
#Simple unit testing
if __name__ == "__main__":
    sdmodel = SDmodel()
    sdmodel.loadSDfile('tester.sdf')
    print sdmodel.length
    sdmodel.setSelected(1)
    print sdmodel.status
    print sdmodel.getMolSvg()[0:100]
    print sdmodel.getMolBlock()[0:100]&lt;/pre&gt;