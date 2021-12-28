# README Bases

- the bases are constructed with the PyMOL Builder Plugin
- protons are added with `h_add` in PyMOL
- hydrogens at O3' and O2P are removed
- protons are renamed manually to conform with the atom names in the AMBER forcefield
- RNA residues are altered from A to RA using PyMOL

- for uridine labeled at C5 (the sugar partial charges are from uridine and the base are from thymidine, the C1' and H1' are then slighly adjusted to reach a charge of -1 for the entire nucleoside)
