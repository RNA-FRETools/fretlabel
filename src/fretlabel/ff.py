#!/usr/bin/env python3

import numpy as np
import pandas as pd
import re
import os
import copy
from biopandas.mol2 import PandasMol2
import json
import pathlib


def pymol_couple_dye2baselinker(
    dye,
    baselinker,
    carbonylC_dye,
    align_baselinker_atoms,
    remove_baselinker_atoms,
    atom_dye="C99",
    atom_baselinker="N99",
):
    """
    Couple a dye to an existing base-linker fragment

    Parameters
    ----------
    dye : str
    baselinker : str
    remove_names : list
        list of atom names to remove
    align_baselinker_atoms : list
        list of atom names of the dye to be used for alignment
    align_baselinker_atoms : list
        list of atom names of the base-linker to be used for alignment
    remove_baselinker_atoms : list
        list of atom names of the base-linker to be removed (i.e. capping group)
    """
    try:
        from pymol import cmd
    except ModuleNotFoundError:
        print("This function can only be executed from within PyMOL")
    else:
        cmd.reinitialize()
        dye = pathlib.Path(dye)
        baselinker = pathlib.Path(baselinker)
        cmd.load(dye)
        cmd.load(baselinker)
        align_dye_atoms = [at.name for at in cmd.get_model(f"neighbor name {carbonylC_dye}", 0).atom]
        align_dye_atoms.sort()
        align_dye_atoms.append(carbonylC_dye)
        cmd.pair_fit(
            *[
                j
                for i in [
                    [
                        "{} and name {}".format(dye.stem, a1),
                        "{} and name {}".format(baselinker.stem, a2),
                    ]
                    for (a1, a2) in zip(align_dye_atoms, align_baselinker_atoms)
                ]
                for j in i
            ]
        )
        cmd.remove("{} and name {}".format(baselinker.stem, "+".join(str(i) for i in remove_baselinker_atoms)))
        cmd.create(f"{dye.stem}_{baselinker.stem}", f"{dye.stem} or {baselinker.stem}")
        cmd.delete(f"{dye.stem} or {baselinker.stem}")
        cmd.bond(f"name {atom_dye}", f"name {atom_baselinker}")
        cmd.alter("all", 'type="ATOM"')
        cmd.set("pdb_use_ter_records", 0)


def pymol_savemol2(filename, molecule, pc_decimals=6, overwrite=False, state=-1):
    """
    Export the selection to a mol2 file with user-defined partial charge precision.
    The PyMOL built-in save function truncates the string to 3 decimals.

    Parameters
    ----------
    pc_decimals : int
        number of decimals for partial charge
    """
    try:
        from pymol import cmd
    except ModuleNotFoundError:
        print("This function can only be executed from within PyMOL")
    else:
        filename = pathlib.Path(filename)
        if filename.is_file() and not overwrite:
            print("File already exists, do not overwrite")
        else:
            # mol2str = cmd.get_str('mol2', molecule, -1)
            model = cmd.get_model(molecule, state)
            substrucutre = []
            visited_res = []
            subst_id = 1
            for i, atom in enumerate(model.atom):
                if atom.resn + str(atom.resi_number) not in visited_res:
                    if atom.flags == 0x08000000:
                        residue_type = "RESIDUE"
                    else:
                        residue_type = "GROUP"
                    substrucutre.append(
                        [
                            subst_id,
                            atom.resn,
                            i + 1,
                            residue_type,
                            atom.resi_number,
                            "****",
                            atom.resn,
                        ]
                    )
                    visited_res.append(atom.resn + str(atom.resi_number))
                    subst_id += 1

            with open(filename, "w") as f:
                f.write("@<TRIPOS>MOLECULE\n{}\n".format(model.molecule.title))
                f.write("{:d} {:d} {:d}\n".format(len(model.atom), len(model.bond), len(visited_res)))
                f.write("SMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n")
                visited_res = []
                subst_id = 0
                for i, atom in enumerate(model.atom):
                    if atom.resn + str(atom.resi_number) not in visited_res:
                        visited_res.append(atom.resn + str(atom.resi_number))
                        subst_id += 1
                    f.write(
                        "{:d}\t{}\t{:0.3f}\t{:0.3f}\t{:0.3f}\t{}\t{:d}\t{}\t{:.{prec}f}\n".format(
                            i + 1,
                            atom.name,
                            *atom.coord,
                            atom.text_type,
                            subst_id,
                            atom.resn,
                            atom.partial_charge,
                            prec=pc_decimals,
                        )
                    )

                f.write("@<TRIPOS>BOND\n")
                for i, bond in enumerate(model.bond):
                    f.write("{:d} {:d} {:d} {:d}\n".format(i + 1, *np.array(bond.index) + 1, bond.order))
                f.write("@<TRIPOS>SUBSTRUCTURE\n")
                for res in substrucutre:
                    f.write("{:d} {} {:d} {} {:d} {} {}\n".format(*res))


def pymol_save_molecule(filename, selection, fmt="mol2", state=-1, overwrite=False):
    """
    Save the molecule

    Parameters
    ----------
    filename : str
    fmt : str
        format of the file (pdb or mol2)
    state : int, optional
        current state = -1
    overwrite : bool
    """
    try:
        from pymol import cmd
    except ModuleNotFoundError:
        print("This function can only be executed from within PyMOL")
    else:
        filename = pathlib.Path(filename)
        if filename.is_file() and not overwrite:
            print("File already exists, do not overwrite")
        else:
            cmd.set("pdb_conect_all", "on")
            cmd.save(filename, selection, state, fmt)


def write_mol2(pandasMol2, filename=None, overwrite=False):
    """
    Write a molecule in the TRIPOS mol2 format

    Parameters
    ----------
    pandasMol2 : biopandas.mol2.pandas_mol2.PandasMol2 instance
    filename : str, optional
    """
    if filename is None:
        filename = "{}_new.mol2".format(pandasMol2.code)
    filename = pathlib.Path(filename)
    atom_str = "@<TRIPOS>ATOM"
    bond_str = "@<TRIPOS>BOND"
    atom_start = pandasMol2.mol2_text.find(atom_str)
    atom_end = atom_start + len(atom_str)
    bond_start = pandasMol2.mol2_text.find(bond_str)

    if filename.is_file() and not overwrite:
        print("File already exists, do not overwrite")
    else:
        with open(filename, "w") as f:
            mol2filestr = "{}\n{}\n{}".format(
                pandasMol2.mol2_text[0:atom_end],
                pandasMol2.df.to_string(header=False, index=False),
                pandasMol2.mol2_text[bond_start:],
            )
            f.write(mol2filestr)


def write_rtp(filename, molecules):
    """
    Write a residue topology parameter (rtp) file

    Parameters
    ----------
    filename : str
    molecule : list of ff.Molecule instances
    """
    if type(molecules) is not list:
        molecules = [molecules]

    with open(filename, "w") as f:
        f.write(
            """[ bondedtypes ]
    ; Col 1: Type of bond
    ; Col 2: Type of angles
    ; Col 3: Type of proper dihedrals
    ; Col 4: Type of improper dihedrals
    ; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.
    ; Col 6: Number of excluded neighbors for nonbonded interactions
    ; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1
    ; Col 8: Remove impropers over the same bond as a proper if it is 1
    ; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih
    1       1          9          4        1         3      1     0\n\n"""
        )
        for mol in molecules:
            f.write("[ {} ]\n".format(mol.moleculetype))
            f.write("[ atoms ]\n")
            f.write(mol.atoms[["atom", "type", "charge", "nr"]].to_string(header=False, index=False))
            f.write("\n[ bonds ]\n")
            f.write(mol.bonds[["ai", "aj"]].to_string(header=False, index=False))
            f.write("\n[ impropers ]\n")
            f.write(mol.impropers[["ai", "aj", "ak", "al"]].to_string(header=False, index=False))
            f.write("\n\n")


def pandasMol2_replace(from_df, to_df, replace_name, subst_value=None, match_name="atom_name"):
    """
    Replace columns from one dataframe to another

    Parameters
    ----------
    from_df : pandas.Dataframe
        DataFrame to take values from
    to_df : pandas.Dataframe
        DataFrame to replace with values
    replace_name : str
        name of the column to be replaced
    subst_value : str
        only replace rows having "subst_value" in the column "subst_name"
    match_name : str
        column to be used as index for dataframe comparison
    """
    if subst_value:
        mask = to_df["subst_name"] == subst_value
    else:
        mask = [True] * to_df.shape[0]
    t1 = to_df.loc[mask, match_name]
    t2 = from_df.set_index(match_name)[replace_name]
    to_df.loc[mask, replace_name] = t1.map(t2)


def update_valency(pandasMol2, bonds_atomNames):
    """
    Update the valency of a set of bonds in a biopandas.mol2.pandas_mol2.PandasMol2 object

    Parameters
    ----------
    pandasMol2 : biopandas.mol2.pandas_mol2.PandasMol2 instance
    bonds_atomNames : pandas.DataFrame
        Dataframe with columns 'name1', 'resn1', 'name2', 'resn2'

    Returns
    -------
    pandasMol2 : biopandas.mol2.pandas_mol2.PandasMol2 instance
    """
    mol2file_multiInd = pandasMol2.df.set_index(["atom_name", "subst_name"])["atom_id"]
    for i, b in bonds_atomNames.iterrows():
        atom_id1 = mol2file_multiInd[b.name1, b.resn1]
        atom_id2 = mol2file_multiInd[b.name2, b.resn2]
        pandasMol2.mol2_text = pandasMol2.mol2_text.replace(
            " {:d} {:d} 1\n".format(atom_id1, atom_id2),
            " {:d} {:d} 2\n".format(atom_id1, atom_id2),
        )
        pandasMol2.mol2_text = pandasMol2.mol2_text.replace(
            " {:d} {:d} 1\n".format(atom_id2, atom_id1),
            " {:d} {:d} 2\n".format(atom_id2, atom_id1),
        )
    return pandasMol2


def check_charge(filename, charge):
    """
    Check the net charge of a mol2 file

    Parameters
    ----------
    filename : str
    charge : float
    """
    mol2 = PandasMol2().read_mol2(str(filename))
    sum_charge = round(mol2._df.charge.sum(), 5)
    if sum_charge == charge:
        print("Check passed!")
    else:
        print("Check failed! The charge is: {:0.4f}".format(sum_charge))


def update_specbond(
    specbond_string,
    inputfile="specbond.dat",
    outputfile="specbond.dat",
    overwrite=False,
):
    """
    Add new special bonds

    Parameters
    ----------
    specbond_string : str
        space-delimited string of following format:
        'resA  atomA  nbondsA  resB  atomB  nbondsB  length  newresA  newresB'
    inputfile : str or path object, optional
    outputfile : str or path object, optional
    overwrite : bool
    """
    inputfile = pathlib.Path(inputfile)
    outputfile = pathlib.Path(outputfile)
    with open(inputfile, "r") as f:
        n_specbonds = int(f.readline())

    specbond_format = [
        "resA",
        "atomA",
        "nbondsA",
        "resB",
        "atomB",
        "nbondsB",
        "length",
        "newresA",
        "newresB",
    ]
    specbonds_df = pd.read_csv(
        inputfile,
        skiprows=1,
        sep="\s+",
        nrows=n_specbonds,
        names=specbond_format,
        na_filter=False,
    )

    try:
        specbond_newline = pd.DataFrame([specbond_string.split()], columns=specbond_format)
    except ValueError:
        print(
            "The specbond string that was passed has a wrong format.\n\nThe format is:\n{}".format(
                " ".join(specbond_format)
            )
        )
    else:
        specbonds_df = specbonds_df.append(specbond_newline).reset_index(drop=True)
        specbonds_df.drop_duplicates(
            inplace=True,
            keep="first",
            subset=["resA", "atomA", "resB", "atomB", "newresA", "newresB"],
        )
        n_specbonds = specbonds_df.shape[0]

        if outputfile.is_file() and not overwrite:
            print("File already exists, do not overwrite. Either change outputfilename or set overwrite=True")
        else:
            with open(outputfile, "w") as f:
                f.write("{:d}\n".format(n_specbonds))
                f.write(specbonds_df.to_string(header=False, index=False) + "\n")


def update_residuetypes(
    residuetypes_string,
    inputfile="residuetypes.dat",
    outputfile="residuetypes.dat",
    overwrite=False,
):
    """
    Add new residue type

    Parameters
    ----------
    residuetypes_string : str
        space-delimited string of following format: 'residue  type'
    inputfile : str or path object, optional
    outputfile : str or path object, optional
    overwrite : bool
    """
    inputfile = pathlib.Path(inputfile)
    outputfile = pathlib.Path(outputfile)
    residuetype_format = ["residue", "type"]
    residuetypes_df = pd.read_csv(inputfile, sep="\s+", names=residuetype_format, na_filter=False)

    try:
        residuetype_newline = pd.DataFrame([residuetypes_string.split()], columns=residuetype_format)
    except ValueError:
        print(
            "The residuetype string that was passed has a wrong format.\n\nThe format is:\n{}".format(
                " ".join(residuetype_format)
            )
        )
    else:
        residuetypes_df = residuetypes_df.append(residuetype_newline).reset_index(drop=True)
        residuetypes_df.drop_duplicates(inplace=True, keep="last", subset=["residue"])

        if outputfile.is_file() and not overwrite:
            print("File already exists, do not overwrite. Either change outputfilename or set overwrite=True")
        else:
            with open(outputfile, "w") as f:
                f.write(residuetypes_df.to_string(header=False, index=False) + "\n")


def update_dye_library(
    dye_entry,
    inputfile="dye_library.json",
    outputfile="dye_library.json",
    overwrite=False,
):
    """
    Add new dye-linker fragment to the FRETlabel library

    Parameters
    ----------
    dye_entry : dict
        dictionary of a dye entry such as:
        "filename":"C3W_DTM", "dye":"sCy3", "base":"DT+RU", "position":"internal"
    inputfile : str or path object, optional
    outputfile : str or path object, optional
    overwrite : bool
    """
    inputfile = pathlib.Path(inputfile)
    outputfile = pathlib.Path(outputfile)
    try:
        with open(inputfile, "r") as f:
            dye_library = json.load(f)
    except FileNotFoundError:
        print("No dye library found. Creating a new one.")
        dye_library = []

    dye_library.append(dye_entry)
    dye_library = pd.DataFrame(dye_library)
    dye_library.drop_duplicates(inplace=True, keep="first")
    dye_library = [val for val in dye_library.T.to_dict().values()]

    if outputfile.is_file() and not overwrite:
        print("File already exists, do not overwrite. Either change outputfilename or set overwrite=True")
    else:
        with open(outputfile, "w") as f:
            json.dump(dye_library, f, indent=2)


class Parameters:
    def __init__(
        self,
        atomtypes,
        bondtypes,
        constrainttypes,
        angletypes,
        propertypes,
        impropertypes,
    ):
        """ """
        self.atomtypes = atomtypes
        self.bondtypes = bondtypes
        self.constrainttypes = constrainttypes
        self.angletypes = angletypes
        self.propertypes = propertypes
        self.impropertypes = impropertypes

    @classmethod
    def read_ff(cls, filelist):
        if isinstance(filelist, str):
            filelist = [filelist]
        fflines = {
            "atomtypes": [],
            "bondtypes": [],
            "constrainttypes": [],
            "angletypes": [],
            "propertypes": [],
            "impropertypes": [],
        }
        for filename in filelist:
            with open(filename, "r") as f:
                lines = f.readlines()

                for i, line in enumerate(lines):
                    match = re.search("\[\s(\w+)\s]", line)
                    if match:
                        key = match.group(1)
                    line = line.lstrip()
                    if line and line[0].isalpha():
                        if (key == "dihedraltypes") and (" 9 " in line):
                            key2 = "propertypes"
                        elif (key == "dihedraltypes") and (" 4 " in line):
                            key2 = "impropertypes"
                        else:
                            key2 = key
                        fflines[key2].append(line)

        atomtypes = pd.DataFrame(
            [x.split() for x in fflines["atomtypes"]],
            columns=["name", "at.num", "mass", "charge", "ptype", "sigma", "epsilon"],
        ).dropna()
        bondtypes = pd.DataFrame(
            [x.split() for x in fflines["bondtypes"]],
            columns=["i", "j", "funct", "b0", "kb"],
        ).dropna()
        constrainttypes = pd.DataFrame(
            [x.split() for x in fflines["constrainttypes"]],
            columns=["i", "j", "funct", "b0"],
        ).dropna()
        angletypes = pd.DataFrame(
            [x.split()[0:6] for x in fflines["angletypes"]],
            columns=["i", "j", "k", "funct", "th0", "cth"],
        ).dropna()
        propertypes = pd.DataFrame(
            [x.split()[0:8] for x in fflines["propertypes"]],
            columns=["i", "j", "k", "l", "funct", "phase", "kb", "pn"],
        ).dropna()
        impropertypes = pd.DataFrame(
            [x.split()[0:8] for x in fflines["impropertypes"]],
            columns=["i", "j", "k", "l", "funct", "phase", "kb", "pn"],
        ).dropna()

        return cls(
            atomtypes,
            bondtypes,
            constrainttypes,
            angletypes,
            propertypes,
            impropertypes,
        )

    @classmethod
    def read_amberdyes(cls, filelist):
        if isinstance(filelist, str):
            filelist = [filelist]
        amberlines = {
            "atomtypes": [],
            "bondtypes": [],
            "constrainttypes": [],
            "angletypes": [],
            "propertypes": [],
            "impropertypes": [],
        }
        for filename in filelist:
            with open(filename, "r") as f:
                lines = f.readlines()

                for i, line in enumerate(lines):
                    match = re.search("\[\s(\w+)\s]", line)
                    if match:
                        key = match.group(1)
                    if "AMBER-DYES" in line:
                        if (key == "dihedraltypes") and (" 9 " in line):
                            key2 = "propertypes"
                        elif (key == "dihedraltypes") and (" 4 " in line):
                            key2 = "impropertypes"
                        else:
                            key2 = key
                        amberlines[key2].append(line)

        atomtypes = (
            pd.DataFrame(
                [x.split() for x in amberlines["atomtypes"]],
                columns=[
                    "name",
                    "at.num",
                    "mass",
                    "charge",
                    "ptype",
                    "sigma",
                    "epsilon",
                    ";",
                    "comment",
                ],
            )
            .dropna()
            .drop(";", axis=1)
        )
        bondtypes = (
            pd.DataFrame(
                [x.split() for x in amberlines["bondtypes"]],
                columns=["i", "j", "funct", "b0", "kb", ";", "comment"],
            )
            .dropna()
            .drop(";", axis=1)
        )
        constrainttypes = (
            pd.DataFrame(
                [x.split() for x in amberlines["constrainttypes"]],
                columns=["i", "j", "funct", "b0", ";", "comment"],
            )
            .dropna()
            .drop(";", axis=1)
        )
        angletypes = (
            pd.DataFrame(
                [x.split() for x in amberlines["angletypes"]],
                columns=["i", "j", "k", "funct", "th0", "cth", ";", "comment"],
            )
            .dropna()
            .drop(";", axis=1)
        )
        propertypes = (
            pd.DataFrame(
                [x.split() for x in amberlines["propertypes"]],
                columns=[
                    "i",
                    "j",
                    "k",
                    "l",
                    "funct",
                    "phase",
                    "kb",
                    "pn",
                    ";",
                    "comment",
                ],
            )
            .dropna()
            .drop(";", axis=1)
        )
        impropertypes = (
            pd.DataFrame(
                [x.split() for x in amberlines["impropertypes"]],
                columns=[
                    "i",
                    "j",
                    "k",
                    "l",
                    "funct",
                    "phase",
                    "kb",
                    "pn",
                    ";",
                    "comment",
                ],
            )
            .dropna()
            .drop(";", axis=1)
        )
        atomtypes["comment"] = "; " + atomtypes["comment"]
        bondtypes["comment"] = "; " + bondtypes["comment"]
        constrainttypes["comment"] = "; " + constrainttypes["comment"]
        angletypes["comment"] = "; " + angletypes["comment"]
        propertypes["comment"] = "; " + propertypes["comment"]
        impropertypes["comment"] = "; " + impropertypes["comment"]
        return cls(
            atomtypes,
            bondtypes,
            constrainttypes,
            angletypes,
            propertypes,
            impropertypes,
        )

    @classmethod
    def read_frcmod(cls, filename, atomtypes_molecule=None):
        """

        # Note: units in AMBER frcmod are different from those in Gromacs!
        # frcmod: A, kcal/mol
        # Gromacs: nm, kJ/mol
        """

        with open(filename, "r") as f:
            lines = f.readlines()
            atomtype_list = []
            bondtype_list = []
            angletype_list = []
            propertype_list = []
            impropertype_list = []

            flag = None
            for i, line in enumerate(lines):
                if "MASS" in line:
                    flag = "atomtype"
                if "BOND" in line:
                    flag = "bondtype"
                if "ANGLE" in line:
                    flag = "angletype"
                if "DIHE" in line:
                    flag = "propertype"
                if "IMPROPER" in line:
                    flag = "impropertype"

                if flag == "atomtype":
                    match = re.search("(\w+\*?)\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
                    if match:
                        atomtype_list.append(match.group(1).strip())
                if flag == "bondtype":
                    match = re.search("(\w+\*?)\s?-(\w+\*?)\s+(\d+\.\d+)\s+(\d+\.\d+)", line)
                    if match:
                        bondtype_list.append(
                            [
                                match.group(1).strip(),
                                match.group(2).strip(),
                                1,
                                float(match.group(4)) / 10,
                                float(match.group(3)) * 4.1868 * 100,
                                "; FRETLABEL",
                            ]
                        )
                if flag == "angletype":
                    match = re.search(
                        "(\w+\*?)\s?-(\w+\*?)\s?-(\w+\*?)\s+(\d+\.\d+)\s+(\d+\.\d+)",
                        line,
                    )
                    if match:
                        angletype_list.append(
                            [
                                match.group(1).strip(),
                                match.group(2).strip(),
                                match.group(3).strip(),
                                1,
                                float(match.group(5)),
                                float(match.group(4)) * 4.1868,
                                "; FRETLABEL",
                            ]
                        )
                if flag == "propertype":
                    match = re.search(
                        "(\w+\s?\*?)-(\w+\s?\*?)-(\w+\s?\*?)-(\w+\*?)\s+(\d)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)",
                        line,
                    )
                    if match:
                        propertype_list.append(
                            [
                                match.group(1).strip(),
                                match.group(2).strip(),
                                match.group(3).strip(),
                                match.group(4).strip(),
                                9,
                                float(match.group(7)),
                                float(match.group(6)) * 4.1868,
                                float(match.group(8)),
                                "; FRETLABEL",
                            ]
                        )
                if flag == "impropertype":
                    match = re.search(
                        "(\w+\s?\*?)-(\w+\s?\*?)-(\w+\s?\*?)-(\w+\*?)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)",
                        line,
                    )
                    if match:
                        impropertype_list.append(
                            [
                                match.group(1).strip(),
                                match.group(2).strip(),
                                match.group(3).strip(),
                                match.group(4).strip(),
                                4,
                                float(match.group(6)),
                                float(match.group(5)) * 4.1868,
                                float(match.group(7)),
                                "; FRETLABEL",
                            ]
                        )

        if (atomtypes_molecule is not None) and atomtype_list:
            atomtypes = atomtypes_molecule[atomtypes_molecule["name"].isin(atomtype_list)]
        else:
            atomtypes = None

        bondtypes = pd.DataFrame(bondtype_list, columns=["i", "j", "funct", "b0", "kb", "comment"])
        angletypes = pd.DataFrame(angletype_list, columns=["i", "j", "k", "funct", "th0", "cth", "comment"])
        propertypes = pd.DataFrame(
            propertype_list,
            columns=["i", "j", "k", "l", "funct", "phase", "kb", "pn", "comment"],
        )
        impropertypes = pd.DataFrame(
            impropertype_list,
            columns=["i", "j", "k", "l", "funct", "phase", "kb", "pn", "comment"],
        )
        if bondtypes.empty:
            bondtypes = None
        if angletypes.empty:
            angletypes = None
        if propertypes.empty:
            propertypes = None
        if impropertypes.empty:
            impropertypes = None

        return cls(atomtypes, bondtypes, None, angletypes, propertypes, impropertypes)

    @classmethod
    def read_specialbond(cls, amberdyes, atoms_amberdyes, atoms_other, ff_name="AMBER-DYES"):
        """
        Alternative constructor of the Parameters class which creates bondtypes, angletypes and propertypes by similarity
        to existing parameters

        Parameters
        ----------
        amberdyes : fretlabel.ff.Parameters instance
        atoms_amberdyes : dict
            dictionary with keys 'bondtypes', 'angletypes' and 'propertypes' and values defining the
        comment : str
            comment string to be added at the end of each line
        """
        bondtype_list = []
        angletype_list = []
        propertype_list = []
        impropertype_list = []
        if "bondtypes" in atoms_amberdyes.keys():
            for a, b in zip(atoms_amberdyes["bondtypes"], atoms_other["bondtypes"]):
                bondtype1 = copy.deepcopy(amberdyes.bondtypes[(amberdyes.bondtypes[["i", "j"]] == a).all(1)])
                bondtype2 = copy.deepcopy(amberdyes.bondtypes[(amberdyes.bondtypes[["i", "j"]] == a[::-1]).all(1)])
                if not bondtype1.empty:
                    bondtype1[["i", "j"]] = b
                    bondtype1["comment"] = "; same as {} {}-{}".format(ff_name, *a)
                    bondtype_list.append(bondtype1)
                if not bondtype2.empty:
                    bondtype2[["i", "j"]] = b[::-1]
                    bondtype2["comment"] = "; same as {} {}-{}".format(ff_name, *a[::-1])
                    bondtype_list.append(bondtype2)
                bondtypes = pd.concat(bondtype_list)
        else:
            bondtypes = None

        if "angletypes" in atoms_amberdyes.keys():
            for a, b in zip(atoms_amberdyes["angletypes"], atoms_other["angletypes"]):
                angletype1 = copy.deepcopy(amberdyes.angletypes[(amberdyes.angletypes[["i", "j", "k"]] == a).all(1)])
                angletype2 = copy.deepcopy(
                    amberdyes.angletypes[(amberdyes.angletypes[["i", "j", "k"]] == a[::-1]).all(1)]
                )
                if not angletype1.empty:
                    angletype1[["i", "j", "k"]] = b
                    angletype1["comment"] = "; same as {} {}-{}-{}".format(ff_name, *a)
                    angletype_list.append(angletype1)
                if not angletype2.empty:
                    angletype2[["i", "j", "k"]] = b[::-1]
                    angletype2["comment"] = "; same as {} {}-{}-{}".format(ff_name, *a[::-1])
                    angletype_list.append(angletype2)
                angletypes = pd.concat(angletype_list)
        else:
            angletypes = None

        if "propertypes" in atoms_amberdyes.keys():
            for a, b in zip(atoms_amberdyes["propertypes"], atoms_other["propertypes"]):
                propertype1 = copy.deepcopy(
                    amberdyes.propertypes[(amberdyes.propertypes[["i", "j", "k", "l"]] == a).all(1)]
                )
                propertype2 = copy.deepcopy(
                    amberdyes.propertypes[(amberdyes.propertypes[["i", "j", "k", "l"]] == a[::-1]).all(1)]
                )
                if not propertype1.empty:
                    propertype1[["i", "j", "k", "l"]] = b
                    propertype1["comment"] = "; same as {} {}-{}-{}-{}".format(ff_name, *a)
                    propertype_list.append(propertype1)
                if not propertype2.empty:
                    propertype2[["i", "j", "k", "l"]] = b[::-1]
                    propertype2["comment"] = "; same as {} {}-{}-{}-{}".format(ff_name, *a[::-1])
                    propertype_list.append(propertype2)
                propertypes = pd.concat(propertype_list)
        else:
            propertypes = None

        if "impropertypes" in atoms_amberdyes.keys():
            for a, b in zip(atoms_amberdyes["impropertypes"], atoms_other["impropertypes"]):
                impropertype1 = copy.deepcopy(
                    amberdyes.impropertypes[(amberdyes.impropertypes[["i", "j", "k", "l"]] == a).all(1)]
                )
                impropertype2 = copy.deepcopy(
                    amberdyes.impropertypes[(amberdyes.impropertypes[["i", "j", "k", "l"]] == a[::-1]).all(1)]
                )
                if not impropertype1.empty:
                    impropertype1[["i", "j", "k", "l"]] = b
                    impropertype1["comment"] = "; same as {} {}-{}-{}-{}".format(ff_name, *a)
                    impropertype_list.append(impropertype1)
                if not impropertype2.empty:
                    impropertype2[["i", "j", "k", "l"]] = b[::-1]
                    impropertype2["comment"] = "; same as {} {}-{}-{}-{}".format(ff_name, *a[::-1])
                    impropertype_list.append(impropertype2)
                impropertypes = pd.concat(impropertype_list)
        else:
            impropertypes = None

        return cls(None, bondtypes, None, angletypes, propertypes, impropertypes)

    def append(self, parameters):
        if self.atomtypes is not None:
            self.atomtypes = self.atomtypes.append(parameters.atomtypes).drop_duplicates()
        if self.bondtypes is not None:
            self.bondtypes = self.bondtypes.append(parameters.bondtypes).drop_duplicates()
        if self.constrainttypes is not None:
            self.constrainttypes = self.constrainttypes.append(parameters.constrainttypes).drop_duplicates()
        if self.angletypes is not None:
            self.angletypes = self.angletypes.append(parameters.angletypes).drop_duplicates()
        if self.propertypes is not None:
            self.propertypes = self.propertypes.append(parameters.propertypes).drop_duplicates()
        if self.impropertypes is not None:
            self.impropertypes = self.impropertypes.append(parameters.impropertypes).drop_duplicates()

    def add2ff(self, ff_folder, outputdir="./"):
        """
        Update the forcefield parameter files ffnonbonded.itp and ffbonded.itp with a new parameter set
        Writes a copy of the forcefield files with the updated parameters to the current working directory

        Parameters
        ----------
        ff_folder : str or path object
            folder where the original ffnonbonded.itp and ffbonded.itp are located
        outputdir :  str or path object
            directory
        """
        ff_folder = pathlib.Path(ff_folder)
        outputdir = pathlib.Path(outputdir)
        for filename in [ff_folder.joinpath("ffnonbonded.itp"), ff_folder.joinpath("ffbonded.itp")]:
            with open(filename, "r") as f:
                lines = f.readlines()
                newlines = ""
                it = enumerate(lines)
                for i, line in it:
                    newlines += line
                    match = re.search("\[\s(\w+)\s]", line)
                    if match:
                        newlines += lines[i + 1]
                        next(it)
                        key = match.group(1)
                        key2 = key
                        if (key == "dihedraltypes") and (" 9 " in lines[i + 2]):
                            key2 = "propertypes"
                        if (key == "dihedraltypes") and (" 4 " in lines[i + 2]):
                            key2 = "impropertypes"

                        df = getattr(self, key2)
                        if df is not None:
                            newlines += df.to_string(header=False, index=False) + "\n"

            with open(outputdir.joinpath(filename.name), "w") as f:
                f.write(newlines)

    def write_atp(self, filename):
        """
        Write a atomtype force field file

        Parameters
        ----------
        filename : str, optional
        """
        with open(filename, "w") as f:
            f.write(self.atomtypes[["name", "mass"]].to_string(header=False, index=False))


class Molecule:
    def __init__(
        self,
        moleculetype,
        atoms,
        bonds,
        angles,
        propers,
        impropers,
        atomtypes_molecule=None,
    ):
        self.moleculetype = moleculetype
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.propers = propers
        self.impropers = impropers
        self.atomtypes_molecule = atomtypes_molecule

    @classmethod
    def read_molecule(cls, filename, comment=None):
        """
        Alternative constructor of the Molecule class which parses an molecule itp file

        Parameters
        ----------
        filename : str
        comment : str
            comment string to be added at the end of each line
        """
        with open(filename, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "[ atomtypes ]" in line:
                    atomtype_ln = i + 1
                if "[ moleculetype ]" in line:
                    moleculetype_ln = i + 1
                    moleculetype = re.search("\w{3}", lines[i + 2]).group()
                if "[ atoms ]" in line:
                    atoms_ln = i + 1
                if "[ bonds ]" in line:
                    bonds_ln = i + 1
                if "[ pairs ]" in line:
                    pairs_ln = i + 1
                if "[ angles ]" in line:
                    angles_ln = i + 1
                if ("[ dihedrals ]" in line) and (" propers" in line):
                    proper_dihedrals_ln = i + 1
                if ("[ dihedrals ]" in line) and (" impropers" in line):
                    improper_dihedrals_ln = i + 1

        atomtypes_molecule = pd.read_csv(
            filename,
            skiprows=atomtype_ln,
            sep="\s+",
            nrows=moleculetype_ln - atomtype_ln - 3,
            comment=";",
            names=["name", "bond_type", "mass", "charge", "ptype", "sigma", "epsilon"],
            na_filter=False,
        )
        atoms = pd.read_csv(
            filename,
            skiprows=atoms_ln,
            sep="\s+",
            nrows=bonds_ln - atoms_ln - 3,
            comment=";",
            names=["nr", "type", "resi", "res", "atom", "cgnr", "charge", "mass"],
            na_filter=False,
        )
        bonds = pd.read_csv(
            filename,
            skiprows=bonds_ln,
            sep="\s+",
            nrows=pairs_ln - bonds_ln - 3,
            comment=";",
            names=["i", "j", "funct", "r", "k"],
            na_filter=False,
        )
        angles = pd.read_csv(
            filename,
            skiprows=angles_ln,
            sep="\s+",
            nrows=proper_dihedrals_ln - angles_ln - 3,
            comment=";",
            names=["i", "j", "k", "funct", "theta", "cth"],
            na_filter=False,
        )
        propers = pd.read_csv(
            filename,
            skiprows=proper_dihedrals_ln,
            sep="\s+",
            nrows=improper_dihedrals_ln - proper_dihedrals_ln - 4,
            comment=";",
            names=["i", "j", "k", "l", "funct", "phase", "kd", "pn"],
            na_filter=False,
        )
        impropers = pd.read_csv(
            filename,
            skiprows=improper_dihedrals_ln,
            sep="\s+",
            comment=";",
            names=["i", "j", "k", "l", "funct", "phase", "kd", "pn"],
            na_filter=False,
        )

        # add masses and atomnumbers to atomtypes_molecule
        masses = {
            "H": 1.00800,
            "C": 12.01000,
            "F": 19.00000,
            "N": 14.01000,
            "O": 16.00000,
            "S": 32.06000,
            "P": 30.97000,
        }
        atomnumbers = {"H": 1, "C": 6, "F": 9, "N": 7, "O": 8, "S": 16, "P": 15}
        for i, row in atomtypes_molecule.iterrows():
            atom = row["name"][0]
            atomtypes_molecule.loc[i, "mass"] = masses[atom]
            atomtypes_molecule.loc[i, "at.num"] = atomnumbers[atom]
        atomtypes_molecule = (
            atomtypes_molecule[
                [
                    "name",
                    "at.num",
                    "bond_type",
                    "mass",
                    "charge",
                    "ptype",
                    "sigma",
                    "epsilon",
                ]
            ]
            .astype({"at.num": int})
            .drop("bond_type", axis=1)
        )

        # map the atom name on the atom number of the bonds, angles and (im)proper dihedrals
        map_corresp = atoms.set_index("nr")["atom"]
        for atm_id in ("i", "j"):
            bonds["a{}".format(atm_id)] = bonds[atm_id].map(map_corresp)
        for atm_id in ("i", "j", "k"):
            angles["a{}".format(atm_id)] = angles[atm_id].map(map_corresp)
        for atm_id in ("i", "j", "k", "l"):
            propers["a{}".format(atm_id)] = propers[atm_id].map(map_corresp)
            impropers["a{}".format(atm_id)] = impropers[atm_id].map(map_corresp)

        if comment is not None:
            bonds["comment"] = "; {}".format(comment)
            angles["comment"] = "; {}".format(comment)
            propers["comment"] = "; {}".format(comment)
            impropers["comment"] = "; {}".format(comment)

        return cls(moleculetype, atoms, bonds, angles, propers, impropers, atomtypes_molecule)

    def change_type(self, name, new_type):
        """
        Change atom type

        Parameters
        ----------
        name : str
        new_type : str
        """
        self.atoms.loc[self.atoms["atom"] == name, "type"] = new_type

    def remove_atom(self, name):
        """
        Remove an atom from a molecule

        Parameters
        ----------
        name : str
        """
        if self.atoms is not None:
            self.atoms = self.atoms.drop(self.atoms[(self.atoms["atom"] == name)].index).reset_index(drop=True)
            self.atoms["nr"] = range(self.atoms.shape[0])
        if self.bonds is not None:
            self.bonds = self.bonds.drop(
                self.bonds[(self.bonds["ai"] == name) | (self.bonds["aj"] == name)].index
            ).reset_index(drop=True)
        if self.angles is not None:
            self.angles = self.angles.drop(
                self.angles[
                    (self.angles["ai"] == name) | (self.angles["aj"] == name) | (self.angles["ak"] == name)
                ].index
            ).reset_index(drop=True)
        if self.propers is not None:
            self.propers = self.propers.drop(
                self.propers[
                    (self.propers["ai"] == name)
                    | (self.propers["aj"] == name)
                    | (self.propers["ak"] == name)
                    | (self.propers["al"] == name)
                ].index
            ).reset_index(drop=True)
        if self.impropers is not None:
            self.impropers = self.impropers.drop(
                self.impropers[
                    (self.impropers["ai"] == name)
                    | (self.impropers["aj"] == name)
                    | (self.impropers["ak"] == name)
                    | (self.impropers["al"] == name)
                ].index
            ).reset_index(drop=True)

    def save_rtp(self, filename):
        write_rtp(filename, self)
