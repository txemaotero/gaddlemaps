"""
This module contains MoleculeTop and AtomTop objects which allows to load
atom and molecule information relative to the bonds between atoms.
"""

import os
from itertools import groupby
from typing import Any, List, Set, Tuple, Generator

from ..parsers import itp_top


class MoleculeTop:
    """
    Loads molecules from a topology file.

    This class behaves like a list of atoms which has bonds defined. The
    appropriate parser will be used based on the input file extension. The
    available parsers are summarized in the class attribute "PARSERS". In this
    attribute, the keys are the files extensions and the values the
    corresponding functions that extracts the information from the files with
    that extensions. These functions should return:

        - The name of the molecule
        - A list with tuples with the atoms and residues names in order of
        appearance in the file.
        - A list with tuples with atoms index (referred to the atoms_info
        indexes) that are bonded.

    Parameters
    ----------
    ftop : string
        The path to the file with the molecule name and bonds information. 
    file_format : str, Optional
        The file extension of ftop. If it is None this will be taken from
        ftop.

    Raises
    ------
    ValueError
        If the file format is not supported.
    IOError
        If the input file misses information.

    """
    PARSERS = {
        '.itp': itp_top,
    }

    def __init__(self, ftop: str, file_format: str = None):
        if file_format is None:
            file_format = os.path.splitext(ftop)[1]
        if file_format not in self.PARSERS:
            raise ValueError(f'{file_format} file format is not supported.')

        self.ftop = ftop
        self.name, atoms_info, atoms_bonds = self.PARSERS[file_format](ftop)
        self.atoms: List['AtomTop'] = []

        for index, atom in enumerate(atoms_info):
            self.atoms.append(AtomTop(*atom, index))
        for bond in atoms_bonds:
            self.atoms[bond[0]].connect(self.atoms[bond[1]])

        if not _are_connected(self.atoms):
            raise IOError(('The molecule is not fully connected.'
                           ' Check your topology file.'))

    def __getitem__(self, index: int) -> 'AtomTop':
        return self.atoms[index]

    def __len__(self) -> int:
        return len(self.atoms)

    def __iter__(self) -> Generator['AtomTop', None, None]:
        for atom in self.atoms:
            yield atom

    def __eq__(self, element: Any) -> bool:
        if isinstance(element, MoleculeTop):
            if element.name == self.name and len(self) == len(element):
                return all(at1 == at2 for at1, at2 in zip(self, element)) # type: ignore
        return False

    def __ne__(self, element: Any) -> bool:
        return not self == element

    def __str__(self) -> str:
        return f'MoleculeTop of {self.name}.'

    __repr__ = __str__

    @property
    def resnames(self) -> List[str]:
        """
        list of str: Residue names of the atoms without consecutive
            repetitions.

        To set this property a list with the same length of residues must be
        passed.
        """
        tot_resnames = ('{:5}{}'.format(atom.resname, atom.resid)
                        for atom in self)  # type: ignore
        return [x[0][:5].strip() for x in groupby(tot_resnames)]
    
    @resnames.setter
    def resnames(self, new_resnames: List[str]):
        if not isinstance(new_resnames, list):
            raise ValueError('new_resnames must be a list of strings.')
        if len(new_resnames) != len(self.resnames):
            raise ValueError((f'Expected {len(self.resnames)} residue name'
                              f' while {len(new_resnames)} given.'))
        resname_index = 0
        actual_resname = self[0].residname
        for atom in self:  # type: ignore
            if atom.residname != actual_resname:
                resname_index += 1
                actual_resname = atom.residname
            atom.resname = new_resnames[resname_index]

    @property
    def resname_len_list(self) -> List[Tuple[str, int]]:
        """
        list of tuple(str, int) : 
            [(resname_1, number_of_atoms_with_resname_1),
            (resname_2, number_of_atoms_with_resname_2), ...]
        """
        tot_resnames = ('{:5}{}'.format(atom.resname, atom.resid)
                        for atom in self)  # type: ignore
        res_len = []
        old_resname: List = []
        for resname in tot_resnames:
            if not old_resname:
                old_resname = [resname, 1]
            elif resname != old_resname[0]:
                res_len.append((old_resname[0][:5].strip(), old_resname[1]))
                old_resname = [resname, 1]
            else:
                old_resname[1] += 1
        res_len.append((old_resname[0][:5].strip(), old_resname[1]))
        return res_len

    def index(self, atom: 'AtomTop') -> int:
        """
        Returns the index of the atom in the molecule.

        Parameters
        ----------
        atom : AtomTop
            The atom to find the index.

        Returns
        -------
        index : int
            The index of the atom in the molecule.

        """
        return self.atoms.index(atom)


class AtomTop:
    """
    Atom with information of its name, residue name and bonds.

    It is also needed the index of the atom in the molecule. The bonds are
    initialized as empty set.

    Parameters
    ----------
    name : str
        Atom name
    resname : str
        Residue name of the atom.
    resid : int
        Residue index of the atom.
    index : int
        Atom index in the molecule.

    Attributes
    ----------
    bonds: set of int
        A set with the hash of the atoms that are connected to self.

    """
    def __init__(self, name: str, resname: str, resid: int, index: int):
        self.name = name
        self.resname = resname
        self.resid = resid
        self.index = index
        self.bonds: Set[int] = set()

    def __repr__(self) -> str:
        return f'Itp atom of {self.name} of molecule {self.resname}'

    __str__ = __repr__

    def __hash__(self) -> int:
        """
        Number in the itp line
        """
        return self.index

    def __eq__(self, atom: Any) -> bool:  # type: ignore
        if isinstance(atom, AtomTop):
            condition = (
                (self.index == atom.index) and
                (self.resname == atom.resname) and
                (self.name == atom.name)
            ) 
            return condition
        return False

    @property
    def residname(self) -> str:
        """
        string: An identifier of the residue (resid+name)
        """
        return '{}{}'.format(self.resid, self.resname)

    def connect(self, atom: 'AtomTop'):
        """
        Connects self with other atom setting the bond.

        Parameters
        ----------
        atom : AtomTop
            Atom to connect.

        """
        self.bonds.add(hash(atom))
        atom.bonds.add(hash(self))

    def closest_atoms(self, natoms: int = 2) -> List[int]:
        """
        Returns a list with natoms index of bonded atoms to self.

        If more than natoms are bonded self, the natoms  with lower id_num are
        returned.

        Parameters
        ----------
        natoms : integer
            The number of atoms to return.

        Returns
        -------
        bonded_atoms : list of int
            The list with the index of natoms atoms bonded self.
        """
        return sorted(self.bonds)[:natoms]


def _are_connected(atoms: List[AtomTop]) -> bool:
    connected_atoms: List[int] = []
    _find_connected_atoms(atoms, 0, connected_atoms)
    return len(connected_atoms) == len(atoms)


def _find_connected_atoms(atoms: List[AtomTop], index: int, connected: list):
     if index not in connected:
         connected.append(index)
     for new_index in atoms[index].bonds:
         if new_index in connected:
             continue
         _find_connected_atoms(atoms, new_index, connected)
