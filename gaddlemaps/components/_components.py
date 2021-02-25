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
'''
This module contains Atom and Molecule objects which connect information from
both gro and itp files.
'''

from collections import defaultdict
from typing import (Any, DefaultDict, Dict, Iterator, List, Optional, Tuple,
                    Union)

import numpy as np
from scipy.spatial.distance import euclidean

from . import AtomGro, AtomTop, MoleculeTop, Residue


class Atom:
    """
    An atom class that wraps the AtomTop and AtomGro classes.

    You can access to the methods and attributes that both AtomGro and
    AtomItp have. To create the atom object, both input atoms should have the
    same resname and name attributes. On the other hand, only the attributes
    from the AtomGro can be changed (e.g. positions, velocities, ...) excluding
    the resname and name.

    Parameters
    ----------
    atom_top : AtomTop
        The AtomTop object.
    atom_gro : AtomGro
        The AtomGro object.

    Raises
    ------
    IOError
        If the atom_gro and atom_top do not correspond to the same atom.
    TypeError
        If the inputs are not instances of the corresponding classes.

    """

    def __init__(self, atom_top: AtomTop, atom_gro: AtomGro):
        if not isinstance(atom_gro, AtomGro):
            raise TypeError('atom_gro input have to be an AtomGro instance.')
        if not isinstance(atom_top, AtomTop):
            raise TypeError('atom_top input have to be an AtomTop instance.')
        if atom_gro.resname != atom_top.resname or atom_gro.name != atom_top.name:
            raise IOError((f'Input atoms do not match:\n-{atom_gro}'
                           f'\n-{atom_top}'))
        self._atom_gro = atom_gro
        self._atom_top = atom_top

    @property
    def atom_gro(self) -> AtomGro:
        """
        AtomGro : The input AtomGro object
        """
        return self._atom_gro

    @property
    def atom_top(self) -> AtomTop:
        """
        AtomTop : The input AtomTop object
        """
        return self._atom_top

    @property
    def top_resid(self) -> int:
        """
        Residue number for the part of the atom with the topology (atom_top).
        """
        return self._atom_top.resid

    @top_resid.setter
    def top_resid(self, new_resid: int):
        self._atom_top.resid = new_resid

    @property
    def gro_resid(self) -> int:
        """
        Residue number for the gro part of the atom (atom_gro).
        """
        return self._atom_gro.resid

    @gro_resid.setter
    def gro_resid(self, new_resid: int):
        self._atom_gro.resid = new_resid

    def __getattr__(self, attr: str) -> Any:
        if attr == 'resid':
            raise AttributeError(('To access resid use gro_resid or'
                                  ' top_resid properties.'))
        if hasattr(self._atom_gro, attr):
            return getattr(self._atom_gro, attr)
        if hasattr(self._atom_top, attr):
            return getattr(self._atom_top, attr)
        raise AttributeError(f'Atom object has no attribute {attr}')

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_atom_top', '_atom_gro']:
            super(Atom, self).__setattr__(attr, value)
        # Special cases that are present in both atoms
        elif attr == 'resid':
            raise AttributeError(('To set resid use gro_resid or'
                                  ' top_resid properties.'))
        elif attr in ('resname', 'name'):
            setattr(self._atom_top, attr, value)
            setattr(self._atom_gro, attr, value)
        elif attr in super(Atom, self).__dir__():
            super(Atom, self).__setattr__(attr, value)
        elif hasattr(self._atom_top, attr):
            setattr(self._atom_top, attr, value)
        elif hasattr(self._atom_gro, attr):
            setattr(self._atom_gro, attr, value)
        else:
            super(Atom, self).__setattr__(attr, value)

    def __str__(self) -> str:
        string = (f'Atom {self.name} of molecule {self.resname}'
                  f' with gro residue number {self.gro_resid}.')
        return string

    __repr__ = __str__

    def __hash__(self) -> int:
        return hash(self._atom_top)

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Atom, self).__dir__())
        _dir.update(dir(self._atom_gro))
        _dir.update(dir(self._atom_top))
        return list(_dir)

    def copy(self) -> 'Atom':
        """
        Returns a copy of self.

        Only the gro atom is copied. The atom_top remains the same.

        Returns
        -------
        new_atom : Atom
            The copied atom.

        """
        return Atom(self._atom_top, self._atom_gro.copy())

    def __eq__(self, atom: Any) -> bool:
        """
        Two atoms are the same if they have the same name, resname, index and top_resid.
        """

        if isinstance(atom, Atom):
            condition = (self.resname == atom.resname and
                         self.name == atom.name and
                         self.index == atom.index and
                         self.top_resid == atom.top_resid)
            return condition
        return False


class Molecule(Residue):
    """
    Loads a molecule combining a MoleculeTop and a list of Residue.

    This class wraps all the features of both MoleculeTop and Residues which
    conform the molecule. When an object is initialized a copy of the input
    residues is stored (to avoid undesired attribute changes). This class
    inherits from Residue so they have the same methods and properties
    (although most of them are reimplemented).

    Parameters
    ----------
    molecule_top : MoleculeTop
        The object with the bonds information of the molecule.
    residues : List of Residue
        An list with the residues that constitute the molecule.

    Raises
    ------
    TypeError
        If the input are instances of wrong type.
    IOError
        If residues do not constitute the molecule_top.

    """
    __excluded__ = ['resname', 'resid', 'residname', 'remove_atom']

    def __init__(self, molecule_top: MoleculeTop, residues: List[Residue]):
        if not _molecule_top_and_residues_match(molecule_top, residues):
            raise IOError(('The molecule can not be initialized. '
                           'The input residues have not the same atoms that '
                           'molecule_top has.'))

        self._molecule_top = molecule_top
        self._residues: List[Residue] = []
        # Save the residue number of each atom
        self._each_atom_resid: List[int] = []

        # Initialize attributes
        for res_index, res in enumerate(residues):
            self._residues.append(res.copy())
            self._each_atom_resid += [res_index] * len(res)

    def __getitem__(self, index: int) -> 'Atom':  # type: ignore
        residue_index = self._each_atom_resid[index]
        atom_index = sum(i == residue_index for i in self._each_atom_resid[:index])
        return Atom(self._molecule_top[index], self._residues[residue_index][atom_index])

    def __iter__(self) -> Iterator['Atom']:  # type: ignore
        index = 0
        for res in self._residues:
            for atom in res:
                yield Atom(self.molecule_top[index], atom)
                index += 1

    def __len__(self) -> int:
        return len(self._each_atom_resid)

    def __str__(self) -> str:
        return f'Molecule of {self.name}.'

    __repr__ = __str__

    def __eq__(self, molecule: Any) -> bool:
        if (isinstance(molecule, Molecule) and
            molecule.name == self.name and
            len(molecule) == len(self)):

            for at1, at2 in zip(self, molecule):
                if at1 != at2:
                    return False
            return True
        return False

    def __ne__(self, molecule: Any) -> bool:
        return not self == molecule

    def __add__(self, other: Union['Residue', 'AtomGro']) -> 'Residue':
        raise NotImplementedError

    def __radd__(self, other: Union['Residue', 'AtomGro']) -> 'Residue':
        raise NotImplementedError

    def __getattribute__(self, attr: str) -> Any:
        """
        Special methods that are implemented for Residue but not for
        Molecule.
        """
        if attr in ['resname', 'resid', 'residname', 'remove_atom']:
            raise AttributeError(attr)
        else:
            return super(Molecule, self).__getattribute__(attr)

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._molecule_top, attr):
            return getattr(self._molecule_top, attr)
        raise AttributeError(f'Molecule object has no attribute {attr}')

    def __setattr__(self, attr: str, value: Any):
        if attr in ['resname', 'resid', 'residname', 'remove_atom']:
            raise AttributeError(attr)
        if attr in ['_molecule_top', '_residues']:
            super(Molecule, self).__setattr__(attr, value)
        elif attr in super(Molecule, self).__dir__():
            super(Molecule, self).__setattr__(attr, value)
        elif hasattr(self._molecule_top, attr):
            setattr(self._molecule_top, attr, value)
        else:
            super(Molecule, self).__setattr__(attr, value)

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Molecule, self).__dir__())
        _dir.update(dir(self._molecule_top))
        return list(_dir)

    @property
    def molecule_top(self) -> MoleculeTop:
        """
        MoleculeTop : The object with the topology information of the molecule.
        """
        return self._molecule_top

    @property
    def residues(self) -> List[Residue]:
        """
        List of Residue : a list with the Residue objects that constitute the
        molecule.
        """
        return self._residues

    @property
    def atoms(self) -> List['Atom']:  # type: ignore
        """
        list of Atom: List with the copies of the atoms of the molecule.
        """
        return [atom.copy() for atom in self]

    @property
    def atoms_positions(self) -> np.ndarray:
        """
        numpy.ndarray((N, 3)) : An array with the atoms positions.
        """
        return np.concatenate([res.atoms_positions for res in self._residues])

    @atoms_positions.setter
    def atoms_positions(self, new_positions: np.ndarray):
        super(Molecule, type(self)).atoms_positions.fset(self, new_positions)

    @property
    def atoms_velocities(self) -> Optional[np.ndarray]:
        """
        numpy.ndarray((N, 3)) or None : An array with the atoms velocities.
            If one of the atoms has no velocity this returns None.
        """
        vels = []
        for res in self._residues:
            res_vel = res.atoms_velocities
            if res_vel is None:
                return None
            vels.append(res_vel)
        return np.concatenate(vels)

    @atoms_velocities.setter
    def atoms_velocities(self, new_velocities: Optional[np.ndarray]):
        super(Molecule, type(self)).atoms_velocities.fset(self, new_velocities) # type: ignore

    @property
    def atoms_ids(self) -> List[int]:
        """
        list of int: A list with the ids of the atoms in the residues.
        """
        return sum((res.atoms_ids for res in self._residues), [])

    @atoms_ids.setter
    def atoms_ids(self, new_ids: List[int]):
        super(Molecule, type(self)).atoms_ids.fset(self, new_ids) # type: ignore

    @property
    def resnames(self) -> List[str]:
        """
        List of string: A list with the names of the residues constituting the
            molecule.

        To set this property, a list of strings with the same length as the
        original must be passed. This will change each residue name. You can
        also pass just a string and this will set all the residue names to the
        same value.
        """
        return [res.resname for res in self._residues]

    @resnames.setter
    def resnames(self, new_resnames: Union[str, List[str]]):
        if isinstance(new_resnames, list) and isinstance(new_resnames[0], str):
            if len(new_resnames) != len(self.resnames):
                raise ValueError(('You should provide a list with'
                                  f' {len(self.resnames)} residues instead of'
                                  f' {len(new_resnames)}.'))

            for atom, res_index in zip(self, self._each_atom_resid):
                atom.resname = new_resnames[res_index]
        elif isinstance(new_resnames, str):
            for atom in self:
                atom.resname = new_resnames
        else:
            raise TypeError('Resnames must be a string or a list of string.')

    @property
    def resids(self) -> List[int]:
        """
        List of int: A list with the residues ids of the residues
            constituting the molecule.

        To set this property, a list of int with the same length as the
        original must be passed. This will change each residue id. You can
        also pass just an int and this will set all the residue ids to the
        same value.
        """
        return [res.resid for res in self._residues]

    @resids.setter
    def resids(self, new_resids: Union[int, List[int]]):
        if isinstance(new_resids, list) and isinstance(new_resids[0], int):
            if len(new_resids) != len(self.resids):
                raise ValueError(('You should provide a list with'
                                  f' {len(self.resids)} residue ids instead'
                                  f' of {len(new_resids)}.'))

            for atom, res_index in zip(self, self._each_atom_resid):
                atom.top_resid = new_resids[res_index]
                atom.gro_resid = new_resids[res_index]
        elif isinstance(new_resids, int):
            for atom in self:
                atom.top_resid = new_resids
                atom.gro_resid = new_resids
        else:
            raise TypeError('Resids must be an int or a list of int.')

    def copy(self, new_residues: List[Residue] = None) -> 'Molecule':
        """
        Returns a copy of the molecule.

        If new_molecule_gro is passed, the old residues will be replaced
        to update the positions. This is used in the extrapolation step.

        NOTE: With this method, the molecule_top used for the Molecule
        initialization remains the same. This means that future changes in
        copied molecules may affect other parts of you code. If you want a
        completely independent new molecule use "deep_copy" method.

        Parameters
        ----------
        new_residues
            List of residues to replace the original positions.

        Returns
        -------
        molecule : Molecule
            The copy of the molecule.

        """
        # Maybe this should be optimized
        if new_residues is None:
            new_residues = self._residues
        return Molecule(self._molecule_top, new_residues)

    def deep_copy(self, new_residues: List[Residue] = None) -> 'Molecule':
        """
        Returns a deep copy of the molecule.

        If new_molecule_gro is passed, the old residues will be replaced
        to update the positions. This is used in the extrapolation step. This
        method generates a new molecule that is not linked to any attribute of
        the original one.

        Parameters
        ----------
        new_residues
            List of residues to replace the original positions.

        Returns
        -------
        molecule : Molecule
            The deep copy of the molecule.

        """
        if new_residues is None:
            new_residues = self._residues
        return Molecule(self._molecule_top.copy(), new_residues)

    def index(self, atom: 'Atom') -> int:
        """
        Returns the index of the atom in the molecule.

        Parameters
        ----------
        atom : Atom
            The atom to find the index.

        Returns
        -------
        index : int
            The index of the atom in the molecule.

        """
        for index, self_atom in enumerate(self):
            if self_atom == atom:
                return index
        raise ValueError(f'{atom} is not in molecule.')

    @property
    def bonds_distance(self) -> Dict[int, List[Tuple[int, float]]]:
        """
        dict of int to list of tuple(int, float): A complex data structure
            that collect the information of the bond distances. The key of the
            property corresponds to the atom index in the molecule. The value
            is a list with tuples. For each tuple, the first value corresponds
            with the index of the bonded atom and the second is the length of
            the bond. This property is used in the alignment process in gaddle
            maps.
        """
        bond_info: DefaultDict[int, List[Tuple[int, float]]] = defaultdict(list)
        for index, atom in enumerate(self):
            for index_to in atom.bonds:
                atom_to = self[index_to]
                distance = euclidean(atom.position, atom_to.position)
                bond_info[index].append((index_to, distance))
        return dict(bond_info)

    @classmethod
    def from_files(cls, fgro: str, ftop: str) -> 'Molecule':
        """
        Loads the molecule from gro and a compatible topology file.

        Parameters
        ----------
        fgro : str
            The file name with the gro.
        ftop : str
            The file name with the top.

        Returns
        -------
        molecule : Molecule
            The initialized molecule.

        """
        from . import System
        syst = System(fgro, ftop)
        if len(syst) == 1:
            return syst[0]
        raise IOError('The input gro file has more than one molecule.')


def _molecule_top_and_residues_match(molecule_top: MoleculeTop,
                                     residues: List[Residue]) -> bool:
    if len(molecule_top) != sum(len(res) for res in residues):
        return False
    index = 0
    for res in residues:
        for atom in res:
            at_top = molecule_top[index]
            if atom.resname != at_top.resname or atom.name != at_top.name:
                return False
            index += 1
    return True
