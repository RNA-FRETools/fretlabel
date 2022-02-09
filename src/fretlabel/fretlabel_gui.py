"""
FRETlabel PyMOL Plugin 

Label nucleic acids with a fluorophores 

(C) Fabio Steffen, University of Zurich
"""

import sys
import pathlib
import json
import numpy as np
import webbrowser

from fretlabel import MODULE_DIR
from fretlabel import __urls__

try:
    from pymol import cmd
    from pymol.Qt import QtWidgets, utils, QtCore
except ModuleNotFoundError:
    print("Pymol is not installed.")

dialog = None


def __init_plugin__(app=None):
    """
    Add FRETlabel plugin to the Plugins Menu
    """
    from pymol.plugins import addmenuitemqt

    addmenuitemqt("FRETlabel", run_plugin_gui)


def run_plugin_gui():
    """
    Create the GUI Window
    """
    global dialog
    if dialog is None:
        dialog = App(_pymol_running=True)
    dialog.show()


class App(QtWidgets.QWidget):
    def __init__(self, _pymol_running=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.uiIcon = str(MODULE_DIR.joinpath("UI", "icon.png"))
        fretlabelUI = str(MODULE_DIR.joinpath("UI", "fretlabel.ui"))
        self.settingsUI = str(MODULE_DIR.joinpath("UI", "settings.ui"))
        self.textUI = str(MODULE_DIR.joinpath("UI", "textpad.ui"))
        utils.loadUi(fretlabelUI, self)
        self.setWindowTitle("FRETlabel")
        self.setWindowIcon(utils.QtGui.QIcon(self.uiIcon))
        self._pymol_running = _pymol_running
        self.docsURL = __urls__["Documentation"]
        self.settingsWindow = QtWidgets.QDialog(self)
        utils.loadUi(self.settingsUI, self.settingsWindow)
        self.settingsWindow.setWindowTitle("FRETlabel - Settings")
        self.settings = {"browser": None, "local_docs": None}
        if MODULE_DIR.joinpath(".fretlabel_settings.conf").is_file():
            with open(MODULE_DIR.joinpath(".fretlabel_settings.conf"), "r") as f:
                self.settings = json.load(f)
        self.textWindow = QtWidgets.QDialog(self)
        self.fileNamePath_pdb = None
        utils.loadUi(self.textUI, self.textWindow)

        # activate / deactivate GUI elements
        self.push_addFragment.setEnabled(False)
        self.push_reloadPDB.setEnabled(False)
        self.push_savePDB.setEnabled(False)
        self.spinBox_atomID.setEnabled(False)

        # add fragments from dye library
        with open(MODULE_DIR.joinpath("dye_library.json"), "r") as f:
            self.dye_lib = json.load(f)
        for frag in self.dye_lib:
            if self.comboBox_selectPosition.findText(frag["position"]) == -1:
                self.comboBox_selectPosition.addItem(frag["position"])
        self.selectChemistry()

        # signals
        self.push_addFragment.clicked.connect(self.addDye)
        self.push_reloadPDB.clicked.connect(self.loadPDBinPyMOL)
        self.push_loadPDB.clicked.connect(self.readPDB)
        self.spinBox_atomID.valueChanged.connect(self.update_atom)
        self.comboBox_chain.currentIndexChanged.connect(self.valid_residues)
        self.push_showText.clicked.connect(self.openPDBFile)
        self.comboBox_selectPosition.currentIndexChanged.connect(self.selectChemistry)
        self.comboBox_selectChemistry.currentIndexChanged.connect(self.selectBase)
        self.comboBox_selectBase.currentIndexChanged.connect(self.selectDye)
        self.comboBox_selectDye.currentIndexChanged.connect(self.valid_residues)
        self.push_demo.clicked.connect(self.runDemo)
        self.push_savePDB.clicked.connect(self.savePDB)
        self.push_docs.clicked.connect(self.openDocumentation)
        self.settingsWindow.push_browser.clicked.connect(self.set_browser)
        self.settingsWindow.push_localdocs.clicked.connect(self.set_localdocsDir)

    def addDye(self):
        """
        Attach the dye to the selected residue
        """
        try:
            cmd.load(
                MODULE_DIR.joinpath(
                    "dyes",
                    "{}.pdb".format(self.fragment["filename"]),
                )
            )
        except cmd.pymol.CmdException:
            print("The selected dye fragment cannot be found")
        else:
            residue_names = self.get_residueNames(self.fragment["filename"])

            # change the residue number of the fragment
            cmd.alter(self.fragment["filename"], "resi={:d}".format(self.resi))
            cmd.alter(self.fragment["filename"], 'chain="{}"'.format(self.chain))
            resin = "chain {} and resn {}+{} and resi \{:d}".format(
                self.chain, self.fragment["base"], self.fragment["linker"], self.resi
            )
            # chain = cmd.get_pdbstr('{} and name C1\''.format(resin))[21]

            # make selections of base and sugar backbone
            nucleic = "resn RA+RG+RC+RU+DA+DG+DC+DT"
            bases = "(name C*+N*+O*+H* and {} and not (name C*'+O*'+O*P*+H*'*+P and {}))".format(resin, resin)
            sugar_backbone = "(name C*'+O*'+O*P*+H*'*+P and {})".format(resin)
            sele_frag_bases = "{} and {} and {}".format(self.fragment["filename"], nucleic, bases)
            sele_pdb_bases = "{} and chain {} and resi \{:d} and {} and {}".format(
                self.fileName_pdb[:-4], self.chain, self.resi, nucleic, bases
            )
            sele_frag_sbb = "{} and {}+{} and {}".format(
                self.fragment["filename"],
                nucleic,
                self.fragment["linker"],
                sugar_backbone,
            )
            sele_pdb_sbb = "{} and chain {} and resi \{:d} and {} and {}".format(
                self.fileName_pdb[:-4], self.chain, self.resi, nucleic, sugar_backbone
            )
            # cmd.select('frag_base', sele_frag_bases)
            # cmd.select('pdb_base', sele_pdb_bases)
            # cmd.select('frag_sbb', sele_frag_sbb)
            # cmd.select('pdb_sbb', sele_pdb_sbb)

            # check NA type of fragment
            if cmd.select("{} and {} and name O2'".format(self.fragment["filename"], nucleic)) > 0:  # RNA
                self.NA_typefrag = "RNA"
            elif cmd.select("{} and {}".format(self.fragment["filename"], nucleic)) > 0:  # DNA
                self.NA_typefrag = "DNA"
            else:
                self.NA_typefrag = None

            if self.fragment["position"] == "internal":
                # (1) align fragment on base of PDB, (2) remove base of PDB,
                # (3) remove sugar-backbone of fragment (4) remove O2' of PDB if NA_typefrag is DNA and alter PDB from RNA to DNA
                # Note: (4) is useful for fragments like DTM which are DNA based but can be inserted at an internal DT or RU
                cmd.align(sele_frag_bases, sele_pdb_bases)
                cmd.remove("{} and {} and {}".format(self.fileName_pdb[:-4], resin, bases))
                cmd.remove("{} and {}".format(self.fragment["filename"], sugar_backbone))
                if self.NA_typefrag == "DNA":
                    cmd.remove(
                        "{} and chain {} and resi \{:d} and name O2'".format(
                            self.fileName_pdb[:-4], self.chain, self.resi
                        )
                    )
                    resn = [r for r in self.fragment["base"].split("+") if "D" in r][0]
                    cmd.alter(
                        "{} and chain {} and resi \{:d}".format(self.fileName_pdb[:-4], self.chain, self.resi),
                        'resn="{}"'.format(resn),
                    )
            else:
                # (1) align fragment on PDB sugar-backbone, (2) extract base from fragment,
                # (3) realign base of fragment on base of PDB, (4) remove entire residue of PDB
                cmd.align(sele_frag_sbb, sele_pdb_sbb)
                cmd.extract("temp_base", sele_frag_bases)
                cmd.align("temp_base", sele_pdb_bases)
                cmd.remove("{} and {}".format(self.fileName_pdb[:-4], resin))
                cmd.create(
                    self.fragment["filename"],
                    "temp_base or {}".format(self.fragment["filename"]),
                )
                cmd.delete("temp_base")

            # add fragment to PDB object
            cmd.create(
                self.fileName_pdb[:-4],
                "{} or {}".format(self.fileName_pdb[:-4], self.fragment["filename"]),
            )
            cmd.delete(self.fragment["filename"])

            # make bond between base and sugar
            if ("A" in self.fragment["base"]) or ("G" in self.fragment["base"]):
                cmd.bond("{} and name N9".format(resin), "{} and name C1'".format(resin))
            else:
                cmd.bond("{} and name N1".format(resin), "{} and name C1'".format(resin))

            # make bond between between consecutive residues at 5'-end/3'-end
            if self.fragment["position"] == "5'-end":
                cmd.bond(
                    "{} and name O3'".format(resin),
                    "chain {} and resi \{:d} and name P".format(self.chain, self.resi + 1),
                )
            elif self.fragment["position"] == "3'-end":
                cmd.bond(
                    "chain {} and resi \{:d} and name O3'".format(self.chain, self.resi - 1),
                    "{} and name P".format(resin),
                )

            self.add_H()

            # rename chain, residues and atoms
            cmd.alter("{} and name OP1".format(resin), 'name="O1P"')
            cmd.alter("{} and name OP2".format(resin), 'name="O2P"')
            linkers = list(np.unique([entry["linker"] for entry in self.dye_lib])) + [
                "DA",
                "DG",
                "DC",
                "DT",
                "RA",
                "RG",
                "RC",
                "RU",
            ]
            for resn in residue_names:
                # cmd.alter('resi \{:d} and resn {}'.format(self.resi, resn), 'chain="{}"'.format(self.chain))
                if resn.strip() in linkers:
                    cmd.alter(
                        "chain {} and resi \{:d} and resn {}".format(self.chain, self.resi, resn),
                        'resn="{}"'.format(self.fragment["filename"][-3:]),
                    )

            # visualization settings
            cmd.show("sticks")
            cmd.color("brightorange", "chain {} and resi \{:d}".format(self.chain, self.resi))
            cmd.zoom(self.fileName_pdb[:-4])

    def add_H(self):
        """
        Add hydrogens to the labeled residue
        """
        if "D" in self.fragment["base"]:
            carbons = ["C1'", "C2'", None, "C3'", "C4'", "C5'", None]
            hydrogens = ["H1'", "H2'1", "H2'2", "H3'", "H4'", "H5'1", "H5'2"]
        else:
            carbons = ["C1'", "C2'", "O2'", "C3'", "C4'", "C5'", None]
            hydrogens = ["H1'", "H2'1", "HO'2", "H3'", "H4'", "H5'1", "H5'2"]
        for c, h in zip(carbons, hydrogens):
            if c is not None:
                cmd.h_add("name {} and resi \{:d} and chain {} and polymer.nucleic".format(c, self.resi, self.chain))
                cmd.alter(
                    "name H01 and resi \{:d} and chain {} and polymer.nucleic".format(self.resi, self.chain),
                    'name="{}"'.format(h),
                )
            else:
                cmd.alter(
                    "name H0* and resi \{:d} and chain {} and polymer.nucleic".format(self.resi, self.chain),
                    'name="{}"'.format(h),
                )

    def loadPDBinPyMOL(self):
        """
        Load a PDB file into PyMOL
        """
        cmd.reinitialize()
        cmd.load(self.fileNamePath_pdb)
        if cmd.select("{} and polymer.nucleic and name O2'".format(self.fileName_pdb[:-4])) > 0:  # RNA
            self.NA_typePDB = "RNA"
        elif cmd.select("{} and polymer.nucleic".format(self.fileName_pdb[:-4])) > 0:  # DNA
            self.NA_typePDB = "DNA"
        else:
            self.NA_typePDB = None
        if cmd.select("{} and resn A+G+C+U".format(self.fileName_pdb[:-4])) > 0:
            self.alter_nucleic(direction="forward")
        else:
            self.pdb_altered = False
        cmd.hide("cartoon")
        cmd.show("sticks")
        self.clean_pdb()
        nucleic_protein = "{} and (polymer.nucleic or polymer.protein)".format(self.fileName_pdb[:-4])
        self.n_atoms = cmd.count_atoms(nucleic_protein)
        self.min_max_residue = self.residue_boundaries(nucleic_protein)
        self.comboBox_chain.clear()
        for chain in cmd.get_chains(nucleic_protein):
            self.comboBox_chain.addItem(chain)
        self.comboBox_chain.setCurrentIndex(0)

    def clean_pdb(self):
        """
        Clean up the PDB file for Gromacs
        """
        cmd.remove("hydrogens and resn A+G+C+U+RA+RG+RC+RU+DA+DG+DC+DU")  # remove hydrogens on nucleic acids
        cmd.remove("inorganic or solvent")
        atom_names = [atom.name for atom in cmd.get_model("resi 1").atom]
        if "O5'" not in atom_names:  # O5' is missing in Rosetta models
            for chain in cmd.get_chains():
                first_resi = int([a.resi for a in cmd.get_model("chain {}".format(chain)).atom][0])
                cmd.edit("chain {} and resi {:d} and name C5'".format(chain, first_resi))
                cmd.attach("O", 2, 2)
                cmd.alter("resi 1 and name O01", 'name="O5\'"')

    def alter_nucleic(self, direction="forward"):
        """
        Alter the nucleic acid residues from A,G,C,U/T into RA,RG,RC,RU (RNA) or DA,DG,DC,DT (DNA) => direction = 'forward'
        Alter the nucleic acid residues from RA,RG,RC,RU (RNA) or DA,DG,DC,DT (DNA) into A,G,C,U/T  => direction = 'backward'
        """
        if direction == "forward":
            if self.NA_typePDB == "RNA":
                for r in ["A", "G", "C", "U"]:
                    cmd.alter(
                        "{} and resn {}".format(self.fileName_pdb[:-4], r),
                        'resn="R{}"'.format(r),
                    )
            else:
                for r in ["A", "G", "C", "T"]:
                    cmd.alter(
                        "{} and resn {}".format(self.fileName_pdb[:-4], r),
                        'resn="D{}"'.format(r),
                    )
            self.pdb_altered = True
        else:
            if self.NA_typePDB == "RNA":
                for r in ["RA", "RG", "RC", "RU"]:
                    cmd.alter(
                        "{} and resn {}".format(self.fileName_pdb[:-4], r),
                        'resn="{}"'.format(r[1]),
                    )
            else:
                for r in ["DA", "DG", "DC", "DU"]:
                    cmd.alter(
                        "{} and resn {}".format(self.fileName_pdb[:-4], r),
                        'resn="{}"'.format(r[1]),
                    )
        cmd.deselect()

    def readPDB(self, fileNamePath_pdb=False):
        """
        Load PDB or CIF file

        Parameters
        ----------
        fileNamePath_pdb : str
        """
        if not fileNamePath_pdb:
            self.fileNamePath_pdb, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Load PDB / CIF", "", "PDB / CIF file (*.pdb *cif);;All Files (*)"
            )
        else:
            self.fileNamePath_pdb = fileNamePath_pdb
        self.fileNamePath_pdb = pathlib.Path(self.fileNamePath_pdb)
        if self.fileNamePath_pdb:
            self.fileName_pdb = self.fileNamePath_pdb.name
            with open(self.fileNamePath_pdb, "r") as f:
                self.pdbText = f.read()
                self.push_showText.setEnabled(True)
            self.loadPDBinPyMOL()
            self.push_reloadPDB.setEnabled(True)
            self.push_savePDB.setEnabled(True)
            self.lineEdit_pdbFile.setText(self.fileName_pdb)

    def selectChemistry(self):
        """
        Update the chemistry dropdown menu
        """
        self.comboBox_selectChemistry.clear()
        currPos = self.comboBox_selectPosition.currentText()
        for frag in self.dye_lib:
            if frag["position"] == currPos:
                if self.comboBox_selectChemistry.findText(frag["chemistry"]) == -1:
                    self.comboBox_selectChemistry.addItem(frag["chemistry"])
        self.selectBase()

    def selectBase(self):
        """
        Update the base dropdown menu
        """
        self.comboBox_selectBase.clear()
        currChem = self.comboBox_selectChemistry.currentText()
        for frag in self.dye_lib:
            if frag["chemistry"] == currChem:
                if self.comboBox_selectBase.findText(frag["base"]) == -1:
                    self.comboBox_selectBase.addItem(frag["base"])
        self.selectDye()

    def selectDye(self):
        """
        Update the dye dropdown menu
        """
        self.comboBox_selectDye.clear()
        currBase = self.comboBox_selectBase.currentText()
        for frag in self.dye_lib:
            if frag["base"] == currBase:
                if self.comboBox_selectDye.findText(frag["dye"]) == -1:
                    self.comboBox_selectDye.addItem(frag["dye"])
        self.valid_residues()

    def valid_residues(self):
        """
        Make list of residues where the selected fragment can be attached to
        """
        currChain = self.comboBox_chain.currentText()
        if currChain:
            currPos = self.comboBox_selectPosition.currentText()
            currChem = self.comboBox_selectChemistry.currentText()
            currBase = self.comboBox_selectBase.currentText()
            currDye = self.comboBox_selectDye.currentText()
            for frag in self.dye_lib:
                if (
                    (frag["position"] == currPos)
                    and (frag["chemistry"] == currChem)
                    and (frag["base"] == currBase)
                    and (frag["dye"] == currDye)
                ):
                    self.fragment = frag
            if self.fileNamePath_pdb:
                if self.fragment["position"] == "5'-end":
                    allowed_resis = " or ".join(
                        [
                            "(chain {} and resi \{})".format(chain, min_max[0])
                            for chain, min_max in self.min_max_residue.items()
                        ]
                    )
                elif self.fragment["position"] == "3'-end":
                    allowed_resis = " or ".join(
                        [
                            "(chain {} and resi \{})".format(chain, min_max[1])
                            for chain, min_max in self.min_max_residue.items()
                        ]
                    )
                else:
                    allowed_resis = " or ".join(
                        [
                            "(chain {} and resi \{}-\{})".format(chain, *min_max)
                            for chain, min_max in self.min_max_residue.items()
                        ]
                    )
                selection = "{} and chain {} and resn {} and ({})".format(
                    self.fileName_pdb[:-4],
                    currChain,
                    self.fragment["base"],
                    allowed_resis,
                )
                pdb_str = cmd.get_pdbstr(selection)
                i = 0
                self.resis = {}
                for chain in self.min_max_residue.keys():
                    self.resis[chain] = []
                    for line in pdb_str.split("\n"):
                        if "CONECT" in line:
                            break
                        try:
                            r = int(line[22:26])
                            if r not in self.resis[chain]:
                                self.resis[chain].append(r)
                            i += 1
                        except (IndexError, ValueError):
                            continue
                if self.resis[currChain]:
                    self.before_resi = min(self.resis[currChain])
                    self.spinBox_atomID.setValue(self.before_resi)
                    self.spinBox_atomID.setMaximum(max(self.resis[currChain]))
                    self.spinBox_atomID.setMinimum(min(self.resis[currChain]))
                    self.spinBox_atomID.setEnabled(True)
                    self.push_addFragment.setEnabled(True)
                else:
                    self.spinBox_atomID.setEnabled(False)
                    self.push_addFragment.setEnabled(False)
                    self.lineEdit_pdbAtom.setText(
                        "no {} at {} in chain {}".format(self.fragment["base"], self.fragment["position"], currChain)
                    )
                self.update_atom()

    def update_atom(self):
        """
        Update the atom edit box and the color of the selected residue in PyMOL
        """
        cmd.color("gray80")
        self.chain = self.comboBox_chain.currentText()
        new_resi = self.spinBox_atomID.value()
        if self.resis[self.chain]:
            if new_resi > self.before_resi:
                while new_resi <= max(self.resis[self.chain]):
                    self.before_resi = new_resi
                    if new_resi in self.resis[self.chain]:
                        self.spinBox_atomID.setValue(new_resi)
                        break
                    else:
                        new_resi += 1
            elif new_resi < self.before_resi:
                while new_resi >= min(self.resis[self.chain]):
                    self.before_resi = new_resi
                    if new_resi in self.resis[self.chain]:
                        self.spinBox_atomID.setValue(new_resi)
                        break
                    else:
                        new_resi -= 1
            self.resi = new_resi
            resn = cmd.get_model("chain {} and resi \{:d}".format(self.chain, new_resi)).atom[0].resn
            resi_str = "{}-{} {:d}".format(self.chain, resn, new_resi)
            self.lineEdit_pdbAtom.setText(resi_str)
            selection = "{} and chain {} and resi \{}".format(self.fileName_pdb[:-4], self.chain, new_resi)
            cmd.color("brightorange", selection)

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
        residue_names = []
        for line in ATOM_str.split("\n"):
            r = line[17:20]
            if r and r not in residue_names:
                residue_names.append(r)
        return residue_names

    def residue_boundaries(self, selection):
        """
        Return the smallest and largest residue number in the selection

        Parameters
        ----------
        selection : str

        Returns
        -------
        min_max_residue : tuple of int
        """
        ATOM_str = cmd.get_pdbstr(selection)
        residue_index = {}
        min_max_residue = {}
        for line in ATOM_str.split("\n"):
            if "CONECT" in line:
                break
            try:
                resi = line[22:26]
                chain = line[21]
            except IndexError:
                continue
            if chain not in residue_index.keys():
                residue_index[chain] = []
            if resi and resi not in residue_index[chain]:
                residue_index[chain].append(int(resi))

        for chain in residue_index.keys():
            min_max_residue[chain] = (
                min(residue_index[chain]),
                max(residue_index[chain]),
            )
        return min_max_residue

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
        filename = MODULE_DIR.joinpath("demo", "DNA.pdb")
        self.readPDB(fileNamePath_pdb=filename)

    def savePDB(self):
        """
        Save the PDB in PDB or CIF format
        """
        if self.pdb_altered == True:
            msg = QtWidgets.QMessageBox()
            if self.NA_typePDB == "RNA":
                message = "Change RNA nucleotides A, G, C and U to RA, RG, RC and RU"
            elif self.NA_typePDB == "DNA":
                message = "Change DNA nucleotides A, G, C and T to RA, RG, RC and RT"
            alter_residues = msg.question(self, "Alter residues", message, msg.Yes | msg.No)
            if alter_residues == msg.No:
                self.alter_nucleic(direction="backward")
        filename, filetype = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save PDB", "", "PDB File (*.pdb);;CIF File (*.cif)"
        )
        if filename:
            cmd.set("pdb_use_ter_records", 0)
            cmd.set("pdb_conect_all", 1)
            cmd.sort(self.fileName_pdb[:-4])
            cmd.save(filename, self.fileName_pdb[:-4], -1, filetype[0:3].lower())
            if self.pdb_altered == True:
                self.alter_nucleic(direction="forward")

    def set_browser(self):
        """
        Define the path to a webbrowser executable (e.g. Firefox, Edge or Chrome)
        """
        browser_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Search for default browser")
        if browser_path:
            self.settingsWindow.lineEdit_browser.setText(browser_path)
            self.settings["browser"] = browser_path
            with open("{}/.fretlabel_settings.conf".format(MODULE_DIR), "w") as f:
                json.dump(self.settings, f, indent=2)

    def set_localdocsDir(self):
        """
        Define path to local docs
        """
        docs_path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select docs directory")
        if docs_path:
            self.settingsWindow.lineEdit_localdocs.setText(docs_path)
            self.settings["local_docs"] = docs_path
            with open("{}/.fretlabel_settings.conf".format(MODULE_DIR), "w") as f:
                json.dump(self.settings, f, indent=2)

    def openDocumentation(self):
        """
        Open documentation in the browser
        """
        try:
            webbrowser.open(self.docsURL)
        except webbrowser.Error:
            if not self.docsURL:
                if not self.settings["local_docs"]:
                    msg = QtWidgets.QMessageBox()
                    msg.setIcon(QtWidgets.QMessageBox.Information)
                    msg.setWindowTitle("Location of docs not configured")
                    msg.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)
                    msg.setText("Press <OK> and specify the path of the local docs.")
                    returnValue = msg.exec_()
                    if returnValue == QtWidgets.QMessageBox.Ok:
                        self.set_localdocsDir()

            if not self.settings["browser"]:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Information)
                msg.setWindowTitle("Web browser not configured")
                msg.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)
                msg.setText(
                    "For the documentation to be displayed, the path to a web browser needs to be configured. Press <OK> and search for the browser executable (Firefox, Chrome or Edge)."
                )
                returnValue = msg.exec_()
                if returnValue == QtWidgets.QMessageBox.Ok:
                    self.set_browser()

            if self.settings["browser"]:
                try:
                    browser = webbrowser.get("{} %s".format(self.settings["browser"]))
                    browser.open("file://{}/index.html".format(self.settings["local_docs"]))
                except webbrowser.Error:
                    self.settings["browser"] = None
                    print("Browser not found!")
                    with open("{}/.fretlabel_settings.conf".format(MODULE_DIR), "w") as f:
                        json.dump(self.settings, f, indent=2)
                except FileNotFoundError:
                    self.settings["local_docs"] = None
                    print("Local docs not found!")
                    with open("{}/.fretlabel_settings.conf".format(MODULE_DIR), "w") as f:
                        json.dump(self.settings, f, indent=2)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = App()
    window.show()
    app.exec_()
