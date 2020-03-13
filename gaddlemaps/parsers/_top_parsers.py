"""
This module defines functions to parse topology files and return a list of the
atoms already bonded.
"""

from typing import Tuple, List, TYPE_CHECKING

from . import ItpLineAtom, ItpFile, ItpLineBonds

if TYPE_CHECKING:
    from ..components import AtomTop


def itp_top(fitp: str) -> Tuple[str, List[Tuple[str, str, int]],
                                List[Tuple[int, int]]]:
    """
    Reads the itp file and returns the information to load a molecule.

    This function extracts the name of the molecule and the atoms and residues
    names. It also extract a list with the bonds.

    Parameters
    ----------
    fitp : str
        The itp file name.

    Returns
    -------
    molecule_name : str
        The name of the molecule
    atoms_info : List of tuples of (str, str, int)
        A list with tuples with the atoms and residues names and residue
        indexes in order of appearance in the file.
    atoms_bonds : List of tuples of (int, int)
        A list with tuples with atoms index (referred to the atoms_info indexes)
        that are bonded.
    
    """
    
    itp_file = ItpFile(fitp)

    if 'moleculetype' not in itp_file:
        raise IOError('The input itp must have "moleculetype" section.')
    if 'atoms' not in itp_file:
        raise IOError('The input itp must have "atoms" section.')

    name = _itp_top_name(itp_file)
    atoms, bonds = _itp_top_atoms(itp_file)
    return name, atoms, bonds


def _itp_top_name(itp_file: ItpFile) -> str:
    for line in itp_file['moleculetype']:
        if line.content:
           return line.name
    raise IOError('There is not molecule name in "moleculetype" sec.')


def _itp_top_atoms(itp_file: ItpFile) -> Tuple[List[Tuple[str, str, int]],
                                               List[Tuple[int, int]]]:
    index = 0
    atoms = []
    atoms_number = {}  # atom number to index
    for atom_line in itp_file['atoms']:
        if atom_line.content:
            atoms.append((atom_line.name, atom_line.resname,
                          atom_line.resid))
            atoms_number[atom_line.number] = index
            index += 1
    if not atoms:
        raise IOError('There are not atoms in the atoms section.')

    # Find bonds in correct indexes
    bonds: List[Tuple[int, int]] = []
    for bond in _parse_itp_bonds(itp_file):
        bonds.append((atoms_number[bond[0]], atoms_number[bond[1]]))
    return atoms, bonds


def _parse_itp_bonds(itp_file: ItpFile) -> List[Tuple[int, int]]:
    condition_constraints = 'constraints' in itp_file
    condition_bonds = 'bonds' in itp_file
    bonds = []
    for key in ('constraints', 'bonds'):
        for bond in itp_file.get(key, []):
            bonds.append((bond.atom_from, bond.atom_to))
    return bonds

