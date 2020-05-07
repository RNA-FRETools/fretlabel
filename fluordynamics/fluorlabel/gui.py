import sys
import os
from pymol import cmd
from pymol.Qt import QtWidgets, utils, QtCore
import json
import pandas as pd
#import webbrowser


#package_directory = os.path.dirname(os.path.abspath(cloud.__file__))

#about = {}
#with open(os.path.join(package_directory, '__about__.py')) as a:
#    exec(a.read(), about)


class App(QtWidgets.QWidget):

    def __init__(self, _pymol_running=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.uiIconPath = '{}/icon'.format(package_directory)
        fluorlabelUI = os.path.join(os.path.dirname(__file__), 'fluorlabel.ui')
        self.textUI = os.path.join(os.path.dirname(__file__), 'textpad.ui')
        utils.loadUi(fluorlabelUI, self)
        self.setWindowTitle("FluorLabel")
        #self.setWindowIcon(utils.QtGui.QIcon(self.uiIconPath))
        self._pymol_running = _pymol_running
        self.readTheDocsURL = None
        self.textWindow = QtWidgets.QDialog(self)
        self.fileNamePath_pdb = None
        utils.loadUi(self.textUI, self.textWindow)
        

        # activate / deactivate GUI elements
        self.push_addFragment.setEnabled(False)
        self.push_reloadPDB.setEnabled(False)
        self.push_savePDB.setEnabled(False)
        self.spinBox_atomID.setEnabled(False)


        # add fragments from dye library
        with open(os.path.join(os.path.dirname(__file__), 'dyes/dye_library.json'), 'r') as f:
            self.dye_lib = json.load(f)
        for frag in self.dye_lib:
            if self.comboBox_selectPosition.findText(frag['position']) == -1:
                self.comboBox_selectPosition.addItem(frag['position'])
        self.selectBase()

        
        # signals
        self.push_addFragment.clicked.connect(self.addDye)
        self.push_reloadPDB.clicked.connect(self.loadPDBinPyMOL)
        self.push_loadPDB.clicked.connect(self.readPDB)
        self.spinBox_atomID.valueChanged.connect(self.update_atom)
        self.push_showText.clicked.connect(self.openPDBFile)
        self.comboBox_selectPosition.currentIndexChanged.connect(self.selectBase)
        self.comboBox_selectBase.currentIndexChanged.connect(self.selectDye)
        self.comboBox_selectDye.currentIndexChanged.connect(self.valid_residues)
        self.push_demo.clicked.connect(self.runDemo)
        self.push_savePDB.clicked.connect(self.savePDB)

    def addDye(self):
        """
        Attach the dye to the selected residue
        """
        try:
            cmd.load(os.path.join(os.path.dirname(__file__), 'dyes/{}.pdb'.format(self.fragment['filename'])))
        except cmd.pymol.CmdException:
            print('The selected dye fragment cannot be found')
        else:
            residue_names = self.get_residueNames(self.fragment['filename'])
            
            # change the residue number of the fragment
            cmd.alter(self.fragment['filename'], "resi={:d}".format(self.resi))
            resin = 'resn {} and resi {:d}'.format(self.fragment['base'], self.resi)
            chain = cmd.get_pdbstr('{} and name C1\''.format(resin))[21]
            
            # make selections of base and sugar backbone
            nucleic = 'resn RA+RG+RC+RU+DA+DG+DC+DT'
            bases = '(name C*+N*+O*+H* and {} and not (name C*\'+O*\'+O*P*+H*\'*+P and {}))'.format(resin,resin)
            sugar_backbone = '(name C*\'+O*\'+O*P*+H*\'*+P and {})'.format(resin)
            sele_frag_bases = '{} and {} and {}'.format(self.fragment['filename'], nucleic, bases)
            sele_pdb_bases = '{} and resi {:d} and {} and {}'.format(self.fileName_pdb[:-4], self.resi, nucleic, bases)
            sele_frag_sbb = '{} and {} and {}'.format(self.fragment['filename'], nucleic, sugar_backbone)
            sele_pdb_sbb = '{} and resi {:d} and {} and {}'.format(self.fileName_pdb[:-4], self.resi, nucleic, sugar_backbone)

            # check NA type of fragment
            if cmd.select('{} and {} and name O2\''.format(self.fragment['filename'], nucleic)) > 0: # RNA
                self.NA_typefrag = 'RNA'
            elif cmd.select('{} and {}'.format(self.fragment['filename'], nucleic)) > 0: # DNA
                self.NA_typefrag = 'DNA'
            else:
                self.NA_typefrag = None

            if self.fragment['position'] == 'internal':
                # (1) align fragment on base of PDB, (2) remove base of PDB, 
                # (3) remove sugar-backbone of fragment (4) remove O2' of PDB if NA_typefrag is DNA and alter PDB from RNA to DNA
                # Note: (4) is useful for fragments like DTM which are DNA based but can be inserted at an internal DT or RU
                cmd.align(sele_frag_bases, sele_pdb_bases)
                cmd.remove('{} and {} and {}'.format(self.fileName_pdb[:-4], resin, bases))
                cmd.remove('{} and {}'.format(self.fragment['filename'], sugar_backbone))
                if self.NA_typefrag == 'DNA':
                    cmd.remove('{} and resi {} and name O2\''.format(self.fileName_pdb[:-4], self.resi))
                    resn = [r for r in self.fragment['base'].split('+') if 'D' in r][0]
                    cmd.alter('{} and resi {}'.format(self.fileName_pdb[:-4], self.resi), 'resn="{}"'.format(resn))
            else:
                # (1) align fragment on PDB sugar-backbone, (2) extract base from fragment, 
                # (3) realign base of fragment on base of PDB, (4) remove entire residue of PDB
                cmd.align(sele_frag_sbb, sele_pdb_sbb)
                cmd.extract('temp_base', sele_frag_bases)
                cmd.align('temp_base', sele_pdb_bases)
                cmd.remove('{} and {}'.format(self.fileName_pdb[:-4], resin))
                cmd.create(self.fragment['filename'], 'temp_base or {}'.format(self.fragment['filename']))
                cmd.remove('temp_base')

            # add fragment to PDB object
            cmd.create(self.fileName_pdb[:-4], '{} or {}'.format(self.fileName_pdb[:-4], self.fragment['filename']))
            cmd.delete(self.fragment['filename'])

            # make bond between base and sugar
            if ('A' in self.fragment['base']) or ('G' in self.fragment['base']):
                cmd.bond('{} and name N9'.format(resin), '{} and name C1\''.format(resin))
            else:
                cmd.bond('{} and name N1'.format(resin), '{} and name C1\''.format(resin))

            # make bond between between consecutive residues at 5'-end/3'end
            if self.fragment['position'] == "5'-end":
                cmd.bond('{} and name O3\''.format(resin), 'resi {} and name P'.format(self.resi+1))
            elif self.fragment['position'] == "3'-end":
                cmd.bond('resi {} and name O3\''.format(self.resi-1), '{} and name P'.format(resin))

            self.add_H()

            # rename chain. residues and atoms
            cmd.alter('{} and name OP1'.format(resin), 'name="O1P"')
            cmd.alter('{} and name OP2'.format(resin), 'name="O2P"')
            print(residue_names)
            for resn in residue_names:
                cmd.alter('resi {:d} and resn {}'.format(self.resi, resn), 'chain="{}"'.format(chain))
                if resn.strip() in ['DA', 'DG', 'DC', 'DT', 'RA', 'RG', 'RC', 'RU', 'POS', 'MLE']:
                    cmd.alter('resi {:d} and resn {}'.format(self.resi, resn), 'resn="{}"'.format(self.fragment['filename'][-3:]))
            
            # visualization settings
            cmd.show('sticks')
            cmd.color('skyblue', 'resi {}'.format(self.resi))       
            cmd.zoom(self.fileName_pdb[:-4])


    def add_H(self):
        """
        Add hydrogens to the labeled residue
        """
        if 'D' in self.fragment['base']:
            carbons = ['C1\'','C2\'',None,'C3\'','C4\'','C5\'',None]
            hydrogens = ['H1\'','H2\'1','H2\'2','H3\'','H4\'','H5\'1','H5\'2']
        else:
            carbons = ['C1\'','C2\'','C3\'','C4\'','C5\'',None]
            hydrogens = ['H1\'','H2','H3\'','H4\'','H5\'1','H5\'2']
        for c,h in zip(carbons,hydrogens):
            if c is not None:
                cmd.h_add('name {} and resi {:d} and polymer.nucleic'.format(c, self.resi))
                cmd.alter('name H01 and resi {:d} and polymer.nucleic'.format(self.resi), 'name="{}"'.format(h))
            else:
                cmd.alter('name H0* and resi {:d} and polymer.nucleic'.format(self.resi), 'name="{}"'.format(h))
        
        
    def loadPDBinPyMOL(self):
        """
        Load a PDB file into PyMOL
        """
        cmd.reinitialize()
        cmd.load(self.fileNamePath_pdb)
        if cmd.select('{} and polymer.nucleic and name O2\''.format(self.fileName_pdb[:-4])) > 0: # RNA
            self.NA_typePDB = 'RNA'
        elif cmd.select('{} and polymer.nucleic'.format(self.fileName_pdb[:-4])) > 0: # DNA
            self.NA_typePDB = 'DNA'
        else:
            self.NA_typePDB = None
        if cmd.select('{} and resn A+G+C+U'.format(self.fileName_pdb[:-4])) > 0:
            self.alter_nucleic(direction='forward')
        else:
            self.pdb_altered = False
        cmd.hide('cartoon')
        cmd.show('sticks')
        nucleic_protein = '{} and not (polymer.nucleic or polymer.protein)'.format(self.fileName_pdb[:-4])
        self.n_atoms = cmd.count_atoms(nucleic_protein)
        self.n_residues = self.count_residues(nucleic_protein)
        self.valid_residues()


    def alter_nucleic(self, direction='forward'):
        """
        Alter the nucleic acid residues from A,G,C,U/T into RA,RG,RC,RU (RNA) or DA,DG,DC,DT (DNA) => direction = 'forward'
        Alter the nucleic acid residues from RA,RG,RC,RU (RNA) or DA,DG,DC,DT (DNA) into A,G,C,U/T  => direction = 'backward'
        """
        if direction == 'forward':
            if self.NA_typePDB == 'RNA':
                for r in ['A', 'G', 'C', 'U']:
                    cmd.alter('{} and resn {}'.format(self.fileName_pdb[:-4], r), 'resn="R{}"'.format(r))
            else:
                for r in ['A', 'G', 'C', 'T']:
                    cmd.alter('{} and resn {}'.format(self.fileName_pdb[:-4], r), 'resn="D{}"'.format(r))
            self.pdb_altered = True
        else:
            if self.NA_typePDB == 'RNA':
                for r in ['RA', 'RG', 'RC', 'RU']:
                    cmd.alter('{} and resn {}'.format(self.fileName_pdb[:-4], r), 'resn="{}"'.format(r[1]))
            else:
                for r in ['DA', 'DG', 'DC', 'DU']:
                    cmd.alter('{} and resn {}'.format(self.fileName_pdb[:-4], r), 'resn="{}"'.format(r[1]))
        cmd.deselect()


    def readPDB(self, fileNamePath_pdb=False):
        """
        Load PDB or CIF file

        Parameters
        ----------
        fileNamePath_pdb : str
        """
        if fileNamePath_pdb is False:
            self.fileNamePath_pdb, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load PDB / CIF', '', "PDB / CIF file (*.pdb *cif);;All Files (*)")
        else:
            self.fileNamePath_pdb = fileNamePath_pdb
        if self.fileNamePath_pdb:
            self.fileName_pdb = self.fileNamePath_pdb.split("/")[-1]
            with open(self.fileNamePath_pdb, 'r') as f:
                self.pdbText = f.read()
                self.push_showText.setEnabled(True)
            self.loadPDBinPyMOL()
            self.push_reloadPDB.setEnabled(True) 
            self.push_savePDB.setEnabled(True)             
            self.lineEdit_pdbFile.setText(self.fileName_pdb)
            

    def selectBase(self):
        """
        Update the base dropdown menu
        """
        self.comboBox_selectBase.clear()
        currPos = self.comboBox_selectPosition.currentText()
        for frag in self.dye_lib:
            if frag['position'] == currPos:
                if self.comboBox_selectBase.findText(frag['base']) == -1:
                    self.comboBox_selectBase.addItem(frag['base'])
        self.selectDye()
                

    def selectDye(self):
        """
        Update the dye dropdown menu
        """
        self.comboBox_selectDye.clear()
        currBase = self.comboBox_selectBase.currentText()
        for frag in self.dye_lib:
            if frag['base'] == currBase:
                if self.comboBox_selectDye.findText(frag['dye']) == -1:
                    self.comboBox_selectDye.addItem(frag['dye'])
        self.valid_residues()


    def valid_residues(self):
        """
        Make list of residues where the selected fragment can be attached to
        """
        currPos = self.comboBox_selectPosition.currentText()
        currBase = self.comboBox_selectBase.currentText()
        currDye = self.comboBox_selectDye.currentText()
        for frag in self.dye_lib:
            if (frag['position']==currPos) and (frag['base']==currBase) and (frag['dye']==currDye):
                self.fragment = frag
        if self.fileNamePath_pdb:
            if self.fragment['position'] == "5'-end":
                allowed_resis = '1'
            elif self.fragment['position'] == "3'-end":
                allowed_resis = str(self.n_residues)
            else:
                allowed_resis = '1-{:d}'.format(self.n_residues)

            selection = '{} and resn {} and resi {}'.format(self.fileName_pdb[:-4], self.fragment['base'], allowed_resis)
            pdb_str = cmd.get_pdbstr(selection)
            i = 0
            self.resis = []
            while i < self.n_atoms:
                try:
                    r = int(pdb_str[22+i*81:26+i*81])
                    if r not in self.resis:
                        self.resis.append(r)
                    i+=1
                except ValueError:
                    break
            if self.resis:
                self.before_resi = min(self.resis)
                self.spinBox_atomID.setValue(self.before_resi)
                self.spinBox_atomID.setMaximum(max(self.resis))
                self.spinBox_atomID.setMinimum(min(self.resis))
                self.spinBox_atomID.setEnabled(True)
                self.push_addFragment.setEnabled(True)
            else:
                self.spinBox_atomID.setEnabled(False)
                self.push_addFragment.setEnabled(False)
                self.lineEdit_pdbAtom.setText('no {} at {}'.format(self.fragment['base'], self.fragment['position']))
            self.update_atom()



    def update_atom(self):
        """
        Update the atom edit box and the color of the selected residue in PyMOL
        """
        cmd.color('gray80')
        new_resi = self.spinBox_atomID.value()
        if self.resis:
            if new_resi > self.before_resi:
                while new_resi <= max(self.resis):
                    self.before_resi=new_resi
                    if new_resi in self.resis:
                        self.spinBox_atomID.setValue(new_resi)
                        break
                    else:
                        new_resi+=1
            elif new_resi < self.before_resi:
                while new_resi >= min(self.resis):
                    self.before_resi=new_resi
                    if new_resi in self.resis:
                        self.spinBox_atomID.setValue(new_resi)
                        break
                    else:
                        new_resi-=1
            self.resi = new_resi
            resi_str = '{}{}'.format(self.fragment['base'], new_resi)
            self.lineEdit_pdbAtom.setText(resi_str)
            selection = '{} and resi {}'.format(self.fileName_pdb[:-4], new_resi)
            cmd.color('skyblue', selection)

    def get_residueNames(self, selection):
        """
        Return a list of all residue names present in the molecule selection

        Parameters
        ----------
        selection : str

        Returns
        -------
        residue_names : list
        """
        ATOM_str = cmd.get_pdbstr(selection)
        residue_names= []
        for line in ATOM_str.split('\n'):
            r = line[17:20]
            if r and r not in residue_names:
                residue_names.append(r)
        return residue_names

    def count_residues(self, selection):
        """
        Count the number of residues in the selection

        Parameters
        ----------
        selection : str

        Returns
        -------
        n_residues : int
        """
        ATOM_str = cmd.get_pdbstr(selection)
        residue_index = []
        for line in ATOM_str.split('\n'):
            resi = line[22:26]
            if resi and resi not in residue_index:
                residue_index.append(resi)
        n_residues = len(residue_index)
        return n_residues


    def openPDBFile(self):
        """
        Show the PDB file as a text file
        """
        self.textWindow.textBrowser_pdbFile.setText(self.pdbText)
        self.textWindow.setWindowTitle("FluorDynamics - {}".format(self.fileName_pdb))
        isOK = self.textWindow.exec_()


    def runDemo(self):
        """
        Display a demo file
        """
        filename = os.path.join(os.path.dirname(__file__), 'demo/p19.pdb')
        self.readPDB(fileNamePath_pdb=filename)

    def savePDB(self):
        """
        Save the PDB in PDB or CIF format
        """
        if self.pdb_altered == True:
            self.alter_nucleic(direction='backward')
        filename, filetype = QtWidgets.QFileDialog.getSaveFileName(self, "Save PDB", "", "PDB File (*.pdb);;CIF File (*.cif)")
        cmd.set('pdb_use_ter_records', 0)
        cmd.sort(self.fileName_pdb[:-4])
        cmd.save(filename, self.fileName_pdb[:-4], -1, filetype[0:3].lower())
        if self.pdb_altered == True:
            self.alter_nucleic(direction='forward')


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = App()
    window.show()
    app.exec_()
