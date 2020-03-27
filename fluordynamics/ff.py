#!/usr/bin/env python3

import numpy as np
import pandas as pd
import re
import os
import copy
from biopandas.mol2 import PandasMol2

from fluordynamics import jupyter

package_directory = os.path.dirname(os.path.relpath(__file__))
fragments_dir=os.path.join(package_directory, 'fragments')

def couple_dye2baselinker2(dye, baselinker, remove_names=None, remove_ids=None, atom_dye='C99', atom_baselinker='N99'):
    """
    Couple a dye to an exisiting base-linker fragment

    Parameters
    ----------
    dye : str
    baselinker : str
    remove_names : list
                   list of atom names to remove
    remove_ids : list
                 list of atom ids to remove
    atom_dye : str
               name of dye atom involved in the bond 
    atom_baselinker : str
                      name of base-linker atom involved in the bond
    """
    cmd_gui = jupyter.connect2pymol()
    try:
        cmd_gui = jupyter.connect2pymol()
    except NameError:
        print('First launch a PyMOL server session with pymol -R')
    else:
        cmd_gui.reinitialize()
        cmd_gui.load('{}/dyes/{}.mol2'.format(fragments_dir, dye))
        cmd_gui.load('{}/base_linkers/{}.mol2'.format(fragments_dir, baselinker))
        if remove_names is not None:
            cmd_gui.remove('{} and name {}'.format(baselinker, '+'.join(str(i) for i in remove_names)))
        if remove_ids is not None:
            cmd_gui.remove('{} and id {}'.format(baselinker, '+'.join(str(i) for i in remove_ids)))
        cmd_gui.fuse('{} and name {}'.format(dye, atom_dye), '{} and name {}'.format(baselinker, atom_baselinker))
        cmd_gui.delete(dye)
        cmd_gui.alter('all', 'type="ATOM"')
        cmd_gui.set('pdb_use_ter_records', 0)
        cmd_gui.set_name(baselinker, '{}_{}'.format(dye, baselinker))

def couple_dye2baselinker(dye, baselinker, align_dye_atoms, align_baselinker_atoms, remove_baselinker_atoms, atom_dye='C99', atom_baselinker='N99'):
    """
    Couple a dye to an exisiting base-linker fragment

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
    cmd_gui = jupyter.connect2pymol()
    try:
        cmd_gui = jupyter.connect2pymol()
    except NameError:
        print('First launch a PyMOL server session with pymol -R')
    else:
        cmd_gui.reinitialize()
        cmd_gui.load('{}/dyes/{}.mol2'.format(fragments_dir, dye))
        cmd_gui.load('{}/base_linkers/{}.mol2'.format(fragments_dir, baselinker))
        cmd_gui.pair_fit(*[j for i in [['{} and name {}'.format(dye, a1), '{} and name {}'.format(baselinker, a2)] for (a1,a2) in zip(align_dye_atoms, align_baselinker_atoms)] for j in i])
        cmd_gui.remove('{} and name {}'.format(baselinker, '+'.join(str(i) for i in remove_baselinker_atoms)))
        cmd_gui.create('{}_{}'.format(dye, baselinker), '{} or {}'.format(dye, baselinker))
        cmd_gui.delete('{} or {}'.format(dye, baselinker))
        cmd_gui.bond('name {}'.format(atom_dye), 'name {}'.format(atom_baselinker))
        cmd_gui.alter('all', 'type="ATOM"')
        cmd_gui.set('pdb_use_ter_records', 0)


def write_mol2(pandasMol2, filename=None, overwrite=False):
    """
    Write a molecule in the TRIPOS mol2 format

    Parameters
    ----------
    pandasMol2 : biopandas.mol2.pandas_mol2.PandasMol2 instance
    filename : str (optional)
    """
    atom_str = '@<TRIPOS>ATOM'
    bond_str = '@<TRIPOS>BOND'
    atom_start = pandasMol2.mol2_text.find(atom_str)
    atom_end = atom_start+len(atom_str)
    bond_start = pandasMol2.mol2_text.find(bond_str)
    
    if filename is None:
        filename = '{}_new.mol2'.format(pandasMol2.code)

    if os.path.isfile(filename) and not overwrite:
        print('File already exists, do not overwrite')
    else:
        with open(filename, 'w') as f:
            mol2filestr = '{}\n{}\n{}'.format(pandasMol2.mol2_text[0:atom_end], pandasMol2.df.to_string(header=False, index=False), pandasMol2.mol2_text[bond_start:])
            f.write(mol2filestr)


def pandasMol2_replace(from_df, to_df, replace_name, subst_value=None, match_name='atom_name'):
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
        mask = to_df['subst_name'] == subst_value
    else:
        mask = [True]*to_df.shape[0]
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
    mol2file_multiInd = pandasMol2.df.set_index(['atom_name','subst_name'])['atom_id']
    for i,b in bonds_atomNames.iterrows():
        atom_id1 = mol2file_multiInd[b.name1, b.resn1]
        atom_id2 = mol2file_multiInd[b.name2, b.resn2]
        pandasMol2.mol2_text = pandasMol2.mol2_text.replace(' {:d} {:d} 1\n'.format(atom_id1, atom_id2), ' {:d} {:d} 2\n'.format(atom_id1, atom_id2))
        pandasMol2.mol2_text = pandasMol2.mol2_text.replace(' {:d} {:d} 1\n'.format(atom_id2, atom_id1), ' {:d} {:d} 2\n'.format(atom_id2, atom_id1))
    return pandasMol2


def update_specbond(specbond_string, filename='specbond.dat', outputdir=None, overwrite=False):
    """
    Add new special bonds

    Parameters
    ----------
    specbond_string : str
                      space-delimited string of following format:
                      'resA  atomA  nbondsA  resB  atomB  nbondsB  length  newresA  newresB'
    filename : str (optional)
    """
    with open(filename, 'r') as f:
        n_specbonds = int(f.readline())

    specbond_format = ['resA', 'atomA', 'nbondsA', 'resB', 'atomB', 'nbondsB', 'length', 'newresA', 'newresB']
    specbonds_df = pd.read_csv(filename, skiprows=1, sep='\s+', nrows=n_specbonds, names=specbond_format, na_filter=False)
    
    try:
        specbond_newline = pd.DataFrame([specbond_string.split()], columns=specbond_format)
    except ValueError:
        print('The specbond string that was passed has a wrong format.\n\nThe format is:\n{}'.format(' '.join(specbond_format)))
    else:
        specbonds_df = specbonds_df.append(specbond_newline).reset_index(drop=True)
        specbonds_df.drop_duplicates(inplace=True, keep='first', subset=['resA', 'atomA', 'resB', 'atomB', 'newresA', 'newresB'])
        n_specbonds = specbonds_df.shape[0]

        
        if outputdir is not None:
            filename = os.path.join(outputdir, filename.split('/')[-1])
        if os.path.isfile(filename):
            filename = filename[:-4]+'_new.dat'
        with open(filename, 'w') as f:
            f.write('{:d}\n'.format(n_specbonds))
            f.write(specbonds_df.to_string(header=False, index=False)+'\n')


def write_ff(ff_folder, amberdyes=None, linker=None, specialbond=None, outputdir='./'):
    """
    Update the forcefield parameter files ffnonbonded.itp and ffbonded.itp with a new parameter set
    Writes a copy of the forcefield files with the updated parameters to the current working directory

    Parameters
    ----------
    ff_folder : str
                folder where the original ffnonbonded.itp and ffbonded.itp are located
    amberdyes : fluordynamics.ff.Parameters instance (optional)
    linker : fluordynamics.ff.Parameters instance (optional)
    specialbond : fluordynamics.ff.Parameters instance (optional)
    """
    for filename in ['{}/ffnonbonded.itp'.format(ff_folder), '{}/ffbonded.itp'.format(ff_folder)]:
        with open(filename, 'r') as f:
            lines = f.readlines()
            keynames = []
            newlines = ''
            it = enumerate(lines)
            for i,line in it:
                newlines += line
                match = re.search('\[\s(\w+)\s]', line)
                if match:
                    newlines += lines[i+1]
                    next(it)
                    key = match.group(1)
                    key2 = key
                    if (key == 'dihedraltypes') and (' 9 ' in lines[i+2]):
                        key2 = 'propertypes'
                    if (key == 'dihedraltypes') and (' 4 ' in lines[i+2]):
                        key2 = 'impropertypes'
                    
                    if amberdyes is not None:
                        amberdyes_df = getattr(amberdyes, key2)
                        if amberdyes_df is not None:
                            newlines += amberdyes_df.to_string(header=False, index=False)+'\n'

                    if linker is not None:
                        linker_df = getattr(linker, key2)
                        if linker_df is not None:
                            newlines += linker_df.to_string(header=False, index=False)+'\n'

                    if specialbond is not None:
                        specialbond_df = getattr(specialbond, key2)
                        if specialbond_df is not None:
                            newlines += specialbond_df.to_string(header=False, index=False)+'\n'                 


        with open(os.path.join(outputdir, filename.split('/')[-1]), 'w') as f:
            f.write(newlines)


class Parameters:
    def __init__(self, atomtypes, bondtypes, constrainttypes, angletypes, propertypes, impropertypes):
        """
        """
        self.atomtypes = atomtypes
        self.bondtypes = bondtypes
        self.constrainttypes = constrainttypes
        self.angletypes = angletypes
        self.propertypes = propertypes
        self.impropertypes = impropertypes


    @classmethod
    def read_amberdyes(cls, filelist):
        if isinstance(filelist, str):
            filelist = [filelist] 
        amberlines = {'atomtypes':[],
                      'bondtypes':[],
                      'constrainttypes':[],
                      'angletypes':[],
                      'propertypes':[],
                      'impropertypes':[]
                      }
        for filename in filelist:
            with open(filename, 'r') as f:
                lines = f.readlines()

                for i,line in enumerate(lines):
                    match = re.search('\[\s(\w+)\s]', line)
                    if match:
                        key = match.group(1)
                    if 'AMBER-DYES' in line:
                        if (key == 'dihedraltypes') and (' 9 ' in line):
                            key2 = 'propertypes'
                        elif (key == 'dihedraltypes') and (' 4 ' in line):
                            key2 = 'impropertypes'
                        else: 
                            key2 = key
                        amberlines[key2].append(line)

        atomtypes = pd.DataFrame([x.split() for x in amberlines['atomtypes']], columns=['name','at.num','mass','charge','ptype','sigma','epsilon', ';','comment']).dropna()
        bondtypes = pd.DataFrame([x.split() for x in amberlines['bondtypes']], columns=['i','j','funct','b0','kb',';', 'comment']).dropna()
        constrainttypes = pd.DataFrame([x.split() for x in amberlines['constrainttypes']], columns=['i','j','funct','b0',';','comment']).dropna()
        angletypes = pd.DataFrame([x.split() for x in amberlines['angletypes']], columns=['i','j','k','funct','th0','cth',';','comment']).dropna()
        propertypes = pd.DataFrame([x.split() for x in amberlines['propertypes']], columns=['i','j','k','l','funct','phase','kd','pn',';','comment']).dropna()
        impropertypes = pd.DataFrame([x.split() for x in amberlines['impropertypes']], columns=['i','j','k','l','funct','phase','kd','pn',';','comment']).dropna()
        return cls(atomtypes, bondtypes, constrainttypes, angletypes, propertypes, impropertypes)


    @classmethod
    def read_frcmod(cls, filename, atomtypes_molecule=None):
        """

        # Note: units in AMBER frcmod are different from those in Gromacs!
        # frcmod: A, kcal/mol
        # Gromacs: nm, kJ/mol
        """

        with open(filename, 'r') as f:
            lines = f.readlines()
            atomtype_list = []
            bondtype_list = []
            angletype_list = []
            propertype_list = []
            impropertype_list = []
            
            flag = None
            for i,line in enumerate(lines):
                if 'MASS' in line:
                    flag = 'atomtype'
                if 'BOND' in line:
                    flag = 'bondtype'   
                if 'ANGLE' in line:
                    flag = 'angletype'  
                if 'DIHE' in line:
                    flag = 'propertype' 
                if 'IMPROPER' in line:
                    flag = 'impropertype' 
                
                if flag == 'atomtype':
                    match = re.search('(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                    if match:
                        atomtype_list.append(match.group(1))
                if flag == 'bondtype':
                    match = re.search('(\w+)\s?-(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                    if match:
                        bondtype_list.append([match.group(1), match.group(2), 1, float(match.group(4))/10, float(match.group(3))*4.1868*100, '; FLUOR-DYNAMICS'])
                if flag == 'angletype':
                    match = re.search('(\w+)\s?-(\w+)\s?-(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                    if match:
                        angletype_list.append([match.group(1), match.group(2), match.group(3), 1, float(match.group(5)), float(match.group(4))*4.1868, '; FLUOR-DYNAMICS'])
                if flag == 'propertype':
                    match = re.search('(\w+\s?\*?)-(\w+\s?\*?)-(\w+\s?\*?)-(\w+\*?)\s+(\d)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                    if match:
                        propertype_list.append([match.group(1), match.group(2), match.group(3), match.group(4), 9, float(match.group(7)), float(match.group(6))*4.1868, float(match.group(8)), '; FLUOR-DYNAMICS'])
                if flag == 'impropertype':
                    match = re.search('(\w+\s?\*?)-(\w+\s?\*?)-(\w+\s?\*?)-(\w+\*?)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                    if match:
                        impropertype_list.append([match.group(1), match.group(2), match.group(3), match.group(4), 4, float(match.group(6)), float(match.group(5))*4.1868, float(match.group(7)), '; FLUOR-DYNAMICS'])

        if (atomtypes_molecule is not None) and atomtype_list:
            atomtypes = atomtypes_molecule[atomtypes_molecule['name'].isin(atomtype_list)]
        else:
            atomtypes = None
        bondtypes = pd.DataFrame(bondtype_list, columns=['i','j','funct', 'b0', 'kb','comment'])
        angletypes = pd.DataFrame(angletype_list, columns=['i','j','k','funct','th0', 'cth','comment'])
        propertypes = pd.DataFrame(propertype_list, columns=['i','j','k','l','funct','phase','kb','pn','comment'])
        impropertypes = pd.DataFrame(impropertype_list, columns=['i','j','k','l','funct','phase','kb','pn','comment'])
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
    def read_specialbond(cls, amberdyes, atoms_amberdyes, atoms_other):
        """
        Alternative constructor of the Parameters class which creates bondtypes, angletypes and propertypes by similarity 
        to existing parameters

        Parameters
        ----------
        amberdyes : fluordynamics.ff.Parameters instance
        atoms_amberdyes : dict
                          dictionary with keys 'bondtypes', 'angletypes' and 'propertypes' 
                          and values defining the 
        comment : str
                  comment string to be added at the end of each line
        """
        bondtype_list = []
        angletype_list = []
        propertype_list = []

        for a,b in zip(atoms_amberdyes['bondtypes'], atoms_other['bondtypes']):
            bondtype1 = copy.deepcopy(amberdyes.bondtypes[(amberdyes.bondtypes[['i','j']] == a).all(1)])
            bondtype2 = copy.deepcopy(amberdyes.bondtypes[(amberdyes.bondtypes[['i','j']] == a[::-1]).all(1)])
            if not bondtype1.empty:
                bondtype1[['i','j']] = b
                bondtype1['comment'] = 'same as AMBER-DYES {}-{}'.format(*a)
                bondtype_list.append(bondtype1)
            if not bondtype2.empty:
                bondtype2[['i','j']] = b[::-1]
                bondtype2['comment'] = 'same as AMBER-DYES {}-{}'.format(*a[::-1])
                bondtype_list.append(bondtype2)
            bondtypes = pd.concat(bondtype_list)
            
                
        for a,b in zip(atoms_amberdyes['angletypes'], atoms_other['angletypes']):
            angletype1 = copy.deepcopy(amberdyes.angletypes[(amberdyes.angletypes[['i','j','k']] == a).all(1)])
            angletype2 = copy.deepcopy(amberdyes.angletypes[(amberdyes.angletypes[['i','j','k']] == a[::-1]).all(1)])
            if not angletype1.empty:
                angletype1[['i','j','k']] = b
                angletype1['comment'] = 'same as AMBER-DYES {}-{}-{}'.format(*a)
                angletype_list.append(angletype1)
            if not angletype2.empty:
                angletype2[['i','j','k']] = b[::-1]
                angletype2['comment'] = 'same as AMBER-DYES {}-{}-{}'.format(*a[::-1])
                angletype_list.append(angletype2)
            angletypes = pd.concat(angletype_list)
                
        for a,b in zip(atoms_amberdyes['propertypes'], atoms_other['propertypes']):
            propertype1 = copy.deepcopy(amberdyes.propertypes[(amberdyes.propertypes[['i','j','k','l']] == a).all(1)])
            propertype2 = copy.deepcopy(amberdyes.propertypes[(amberdyes.propertypes[['i','j','k','l']] == a[::-1]).all(1)])
            if not propertype1.empty:
                propertype1[['i','j','k','l']] = b
                propertype1['comment'] = 'same as AMBER-DYES {}-{}-{}-{}'.format(*a)
                propertype_list.append(propertype1)
            if not propertype2.empty:
                propertype2[['i','j','k','l']] = b[::-1]
                propertype2['comment'] = 'same as AMBER-DYES {}-{}-{}-{}'.format(*a[::-1])
                propertype_list.append(propertype2)
            propertypes = pd.concat(propertype_list)
        return cls(None, bondtypes, None, angletypes, propertypes, None)


    def write_atp(self, filename):
        """
        Write a atomtype force field file

        Parameters
        ----------
        filename : str (optional)
        """
        with open(filename, 'w') as f:
            f.write(self.atomtypes[['name', 'mass']].to_string(header=False, index=False))


class Molecule:
    def __init__(self, moleculetype, atoms, bonds, angles, propers, impropers, atomtypes_molecule=None):
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
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i,line in enumerate(lines):
                if '[ atomtypes ]' in line:
                    atomtype_ln = i+1
                if '[ moleculetype ]' in line:
                    moleculetype_ln = i+1
                    moleculetype = re.search('\w{3}',lines[i+2]).group()
                if '[ atoms ]' in line:
                    atoms_ln = i+1
                if '[ bonds ]' in line:
                    bonds_ln = i+1
                if '[ pairs ]' in line:
                    pairs_ln = i+1
                if '[ angles ]' in line:
                    angles_ln = i+1  
                if ('[ dihedrals ]' in line) and (' propers' in line):
                    proper_dihedrals_ln = i+1 
                if ('[ dihedrals ]' in line) and (' impropers' in line):
                    improper_dihedrals_ln = i+1 

        atomtypes_molecule = pd.read_csv(filename, skiprows=atomtype_ln, sep='\s+', nrows=moleculetype_ln-atomtype_ln-3, comment=';', names=['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon'], na_filter=False)     
        atoms = pd.read_csv(filename, skiprows=atoms_ln, sep='\s+', nrows=bonds_ln-atoms_ln-3, comment=';', names=['nr', 'type', 'resi', 'res', 'atom', 'cgnr', 'charge', 'mass'], na_filter=False)
        bonds = pd.read_csv(filename, skiprows=bonds_ln, sep='\s+', nrows=pairs_ln-bonds_ln-3, comment=';', names=['i', 'j', 'funct', 'r', 'k'], na_filter=False)
        angles = pd.read_csv(filename, skiprows=angles_ln, sep='\s+', nrows=proper_dihedrals_ln-angles_ln-3, comment=';', names=['i', 'j', 'k', 'funct', 'theta', 'cth'], na_filter=False)
        propers = pd.read_csv(filename, skiprows=proper_dihedrals_ln, sep='\s+', nrows=improper_dihedrals_ln-proper_dihedrals_ln-4, comment=';', names=['i', 'j', 'k', 'l', 'funct', 'phase', 'kd', 'pn'], na_filter=False)
        impropers = pd.read_csv(filename, skiprows=improper_dihedrals_ln, sep='\s+', comment=';', names=['i', 'j', 'k', 'l', 'funct', 'phase', 'kd', 'pn'], na_filter=False)

        # add masses and atomnumbers to atomtypes_molecule
        masses = {'H':1.00800, 'C': 12.01000, 'F':19.00000, 'N':14.01000, 
                   'O':16.00000, 'S':32.06000, 'P':30.97000}
        atomnumbers = {'H':1, 'C':6 , 'F':9, 'N':7, 
                        'O':8, 'S':16, 'P':15}
        for i,row in atomtypes_molecule.iterrows():
            atom = row['name'][0]
            atomtypes_molecule.loc[i,'mass'] = masses[atom]
            atomtypes_molecule.loc[i,'at.num'] = atomnumbers[atom]
        atomtypes_molecule = atomtypes_molecule[['name', 'at.num', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']].astype({'at.num': int}).drop('bond_type', axis=1)

        # map the atom name on the atom number of the bonds, angles and (im)proper dihedrals
        map_corresp = atoms.set_index('nr')['atom']
        for atm_id in ('i','j'):
            bonds[atm_id] = bonds[atm_id].map(map_corresp)
        for atm_id in ('i','j','k'):
            angles[atm_id] = angles[atm_id].map(map_corresp)
        for atm_id in ('i','j','k','l'):
            propers[atm_id] = propers[atm_id].map(map_corresp)
            impropers[atm_id] = impropers[atm_id].map(map_corresp)

        if comment is not None:
            bonds['comment'] = '; {}'.format(comment)
            angles['comment'] = '; {}'.format(comment)
            propers['comment'] = '; {}'.format(comment)
            impropers['comment'] = '; {}'.format(comment)

        return cls(moleculetype, atoms, bonds, angles, propers, impropers, atomtypes_molecule)

    def change_type(self, name, new_type):
        """
        Change atom type

        Parameters
        ----------
        name : str
        new_type : str
        """
        self.atoms.loc[self.atoms['atom'] == name, 'type'] = new_type

    def remove_atom(self, name):
        """
        Remove an atom from a molecule

        Parameters
        ----------
        name : str
        """
        self.atoms = self.atoms.drop(self.atoms[self.atoms['atom'] == name].index).reset_index()
        self.atoms['nr'] = range(self.atoms.shape[0])


    def write_rtp(self, filename):
        """
        Write a molecule in the TRIPOS mol2 format

        Parameters
        ----------
        filename : str (optional)
        """
        with open(filename, 'w') as f:
            f.write("""[ bondedtypes ]
        ; Col 1: Type of bond
        ; Col 2: Type of angles
        ; Col 3: Type of proper dihedrals
        ; Col 4: Type of improper dihedrals
        ; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.
        ; Col 6: Number of excluded neighbors for nonbonded interactions
        ; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1
        ; Col 8: Remove impropers over the same bond as a proper if it is 1
        ; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih
        1       1          9          4        1         3      1     0\n\n""")
            f.write('[ {} ]\n'.format(self.moleculetype))
            f.write('[ atoms ]\n')
            f.write(self.atoms[['atom', 'type', 'charge', 'nr']].to_string(header=False, index=False))
            f.write('\n[ bonds ]\n')
            f.write(self.bonds[['i', 'j']].to_string(header=False, index=False))
            f.write('\n[ impropers ]\n')
            f.write(self.impropers[['i', 'j', 'k', 'l']].to_string(header=False, index=False))
            f.write('\n')


