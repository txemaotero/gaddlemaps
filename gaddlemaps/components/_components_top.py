#    Gaddlemaps python module.
#    Copyright (C) 2019-2021 José Manuel Otero Mato, Hadrián Montes Campos, Luis Miguel Varela Cabo
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
This module contains MoleculeTop and AtomTop objects which allows to load
atom and molecule information relative to the bonds between atoms.
"""

import os
from itertools import groupby
from typing import Any, Generator, List, Set, Tuple, overload

from ..parsers import read_topology


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

    def __init__(self, ftop: str, file_format: str = None):

        self.ftop = ftop
        self.name, atoms_info, atoms_bonds = read_topology(ftop,
                                                           file_format=file_format)
        self.atoms: List['AtomTop'] = []

        for index, atom in enumerate(atoms_info):
            self.atoms.append(AtomTop(*atom, index))
        for bond in atoms_bonds:
            self.atoms[bond[0]].connect(self.atoms[bond[1]])

    @overload
    def __getitem__(self, index: int) -> 'AtomTop':
        ...
    @overload
    def __getitem__(self, index: slice) -> List['AtomTop']:
        ...
    def __getitem__(self, index):
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
        for atom in self:
            if atom.residname != actual_resname:
                resname_index += 1
                actual_resname = atom.residname
            atom.resname = new_resnames[resname_index]

    @property
    def resids(self) -> List[int]:
        """
        list of str: Residue names of the atoms without consecutive
            repetitions.

        To set this property a list with the same length of residues must be
        passed.
        """
        tot_resids = ('{:5}{}'.format(atom.resname, atom.resid)
                        for atom in self)  # type: ignore
        return [int(x[0][5:].strip()) for x in groupby(tot_resids)]

    @resids.setter
    def resids(self, new_resids: List[int]):
        if not isinstance(new_resids, list):
            raise ValueError('new_resids must be a list of strings.')
        if len(new_resids) != len(self.resids):
            raise ValueError((f'Expected {len(self.resids)} residue name'
                              f' while {len(new_resids)} given.'))
        resid_index = 0
        actual_resid = self[0].residname
        for atom in self:  # type: ignore
            if atom.residname != actual_resid:
                resid_index += 1
                actual_resid = atom.residname
            atom.resid = new_resids[resid_index]

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

    def copy(self) -> 'MoleculeTop':
        """
        Returns a copy of the molecule_top.

        The atoms forming the copy are not the same objects as the original
        molecule so you do not have to worry about linked objects.

        Returns
        -------
        molecule_top : MoleculeTop
            The copy of the molecule.
        """
        mol = self.__new__(self.__class__)
        mol.ftop = self.ftop
        mol.name = self.name
        mol.atoms = [atom.copy() for atom in self]
        return mol


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
                (self.name == atom.name) and
                (self.bonds == atom.bonds)
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

    def copy(self) -> 'AtomTop':
        """
        Returns a copy of the current atom.

        Returns
        -------
        atom_top : AtomTop
            The copy of the atom.
        """
        atom = AtomTop(self.name, self.resname, self.resid, self.index)
        atom.bonds = self.bonds.copy()
        return atom
