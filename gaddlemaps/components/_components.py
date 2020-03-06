# -*- coding: utf-8 -*-
'''
This module contains Atom and Molecule objects which connect information from
both gro and itp files.
'''

import warnings
import os
from collections import defaultdict
from typing import Optional, Any, List, Dict, Tuple, Union, Iterator
from scipy.spatial.distance import euclidean
from . import (AtomGro, AtomItp, GeneralAtom, GeneralMolecule, MoleculeItp,
               MoleculeGro, MacromoleculeGro, Residue)
from ..parsers import GroLine

class InfoDict():
    def __init__(self):
        self.molecule_itp: str = ""
        self.molecule_gro: List[Groline] = []
    
    def to_dict(self) -> Dict[str, Any]:
        return {"molecule_itp": self.molecule_itp,
                "molecule_gro": self.molecule_gro}
        
    def __dict__(self):
        return self.to_dict()
    
    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> 'InfoDict':
        out = InfoDict()
        
        molecule_itp = value["molecule_itp"]
        if not isinstance(molecule_itp, str):
            raise ValueError("The value molecule_itp must be a string")
        
        out.molecule_itp = molecule_itp
        
        out.molecule_gro = value["molecule_gro"]  # type: ignore
        return out
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, InfoDict):
            return False
        return self.to_dict() == other.to_dict()
        


class Atom:
    """
    An atom class that wraps the AtomItp and AtomGro classes.

    You can access to the methods and attributes that both AtomGro and
    AtomItp have. To create the atom object, both input atoms should have the
    same resname and atomname attributes. On the other hand, only the attributes
    from the AtomGro can be changed (e.g. positions, velocities, ...) excluding
    the resname and atomname. 

    Parameters
    ----------
    atom_itp : AtomItp
        The AtomItp object.
    atom_gro : AtomGro
        The AtomGro object.

    Raises
    ------
    IOError
        If the atom_gro and atom_itp do not correspond to the same atom.
    TypeError
        If the inputs are not instances of the corresponding classes.

    """

    def __init__(self, atom_itp: AtomItp, atom_gro: AtomGro):
        super(Atom, self).__init__()
        if not isinstance(atom_gro, AtomGro):
            raise TypeError('atom_gro input have to be an AtomGro instance.')
        if not isinstance(atom_itp, AtomItp):
            raise TypeError('atom_itp input have to be an AtomItp instance.')
        if atom_gro.residname != atom_itp.residname:
            raise IOError((f'Input atoms do not match:\n-{atom_gro}'
                           '\n-{atom_itp}'))
        self._atom_gro = atom_gro
        self._atom_itp = atom_itp

    @property
    def atom_gro(self) -> AtomGro:
        """
        AtomGro : The input AtomGro object
        """
        return self._atom_gro

    @property
    def atom_itp(self) -> AtomItp:
        """
        AtomItp : The input AtomItp object
        """
        return self._atom_itp

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._atom_gro, attr):
            return getattr(self._atom_gro, attr)
        elif hasattr(self._atom_itp, attr):
            return getattr(self._atom_itp, attr)
        else:
            raise AttributeError(f'Atom object has no attribute {attr}')

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_atom_itp', '_atom_gro']:
            super(Atom, self).__setattr__(attr, value)
        # Special cases that are present in both atoms
        elif attr in ('resname', 'atomname'):
            raise AttributeError(f'Attribute {attr} can not be changed')
        elif attr in super(Atom, self).__dir__():
            super(Atom, self).__setattr__(attr, value)
        elif hasattr(self._atom_itp, attr):
            setattr(self._atom_itp, attr, value)
        elif hasattr(self._atom_gro, attr):
            setattr(self._atom_gro, attr, value)
        else:
            super(Atom, self).__setattr__(attr, value)

    def __str__(self) -> str:
        string = (f'Atom {self.atomname} of molecule {self.resname}'
                  f' with residue number {self.resid}.')
        return string

    __repr__ = __str__

    def __hash__(self) -> int:
        return hash(self._atom_itp)

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Atom, self).__dir__())
        _dir.update(dir(self._atom_gro))
        _dir.update(dir(self._atom_itp))
        return list(_dir)

    def copy(self) -> 'Atom':
        """
        Returns a copy of self.

        Only the gro atom is copied. The itp one remains the same.

        Returns
        -------
        new_atom : Atom
            The copied atom.

        """
        return Atom(self._atom_gro.copy(), self._atom_itp)

    def __eq__(self, atom: Any) -> bool:
        if isinstance(atom, Atom):
            return self.residname == atom.residname
        return False


class Molecule:
    """
    Loads a molecule combining a MoleculeItp and a list of Residue.

    This class wraps all the features of both MoleculeItp and Residues which
    conforms the molecule. When an object is initialized a copy of the input
    residues is stored (to avoid undesired attribute changes).

    Parameters
    ----------
    molecule_itp : MoleculeItp
        The object with the .itp information of the molecule.
    residues :  of Residue
        An iterable (e.g. list or tuple) with the residues that constitute the
        molecule.

    Raises
    ------
    TypeError
        If the input are instances of wrong type.
    IOError
        If residues do not constitute the molecule_itp.

    """

    def __init__(self, molecule_itp: MoleculeItp, residues: List[Residue]):
        if not _molecule_itp_and_residues_match(molecule_itp, residues):
            raise IOError(('The molecule can not be initialized. '
                           'The input residues have not the same atoms that '
                           'molecule_itp has.'))

        self._molecule_itp = molecule_itp
        self._residues: List[Residue] = []
        self._atoms: List[Atom] = []
        # Save the residue number of each atom
        self._each_atom_resid: List[int] = []

        # Initialize attributes
        index = 0
        for res_index, res in enumerate(residues):
            self._residues.append(res.copy())
            self._each_atom_resid.append(res_index)
            for atom in res:  # type: ignore
                self._atoms.append(Atom(molecule_itp[index], atom))
                index += 1

    @property
    def molecule_itp(self) -> MoleculeItp:
        """
        MoleculeItp : The object with the .itp information of the molecule.
        """
        return self._molecule_itp

    @property
    def residues(self) -> List[Residue]:
        """
        List of Residue : a list with the Residue objects that constitute the
        molecule.
        """
        return self._residues

    def __getitem__(self, index: int) -> 'Atom':
        return self._atoms[index]

    def __len__(self) -> int:
        return len(self._atoms)

    def __iter__(self) -> Iterator['Atom']:
        for atom in self._atoms:
            yield atom

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._molecule_itp, attr):
            return getattr(self._molecule_itp, attr)
        raise AttributeError(f'Molecule object has no attribute {attr}')

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_molecule_itp', '_residues']:
            super(Molecule, self).__setattr__(attr, value)
        elif attr in super(Molecule, self).__dir__():
            super(Molecule, self).__setattr__(attr, value)
        elif hasattr(self._molecule_itp, attr):
            setattr(self._molecule_itp, attr, value)
        else:
            super(Molecule, self).__setattr__(attr, value)

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

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Molecule, self).__dir__())
        _dir.update(dir(self._molecule_itp))
        return list(_dir)

    @property
    def atoms(self) -> List['Atom']:
        """
        list of Atom: List with the copies of the atoms of the molecule.
        """
        return [atom.copy() for atom in self._atoms]

    @property
    def resnames(self) -> List[str]:
        """
        List of string: A list with the names of the residues constituting the
            molecule.
        """
        return [res.resname for res in self._residues]

    @resname.setter
    def resnames(self, new_resnames: List[str]):
        if isinstance(new_resnames, list) and isinstance(new_resnames, str):
            if len(new_resnames) != len(self.resnames):
                raise ValueError(('You should provide a list with'
                                  f' {len(self.resnames)} residues instead of'
                                  f' {len(new_resnames)}.'))

            for atom, res_index in zip(self, self._each_atom_resid):
                atom.resname = new_resnames[res_index]
            for res, resname in zip(self._residues, new_resnames):
                res.resname = resname
        else:
            raise TypeError('Resname must be a string.')
    # TODO: Aquí me quedo
    def copy(self, new_molecule_gro: Optional[MoleculeGro] = None) -> 'Molecule':
        """
        Returns a copy of the molecule.

        If new_molecule_gro is passed, the old molecule_gro will be replaced
        to update the positions.

        Parameters
        ----------
        new_molecule_gro: MoleculeGro
            Molecule to replace the original positions.

        Returns
        -------
        molecule : MoleculeGro
            The copy of the molecule.

        """
        mol = self.__new__(self.__class__)
        mol._molecule_itp = self._molecule_itp
        if new_molecule_gro is None:
            mol._molecule_gro = self.molecule_gro.copy()
        else:
            if GeneralMolecule.are_the_same_molecule(self.molecule_gro,
                                                     new_molecule_gro):
                mol._molecule_gro = new_molecule_gro.copy()
            else:
                raise IOError(('The input molecule_gro can not replace the old'
                               ' one, please try to create a new instance.'))
        mol._atoms = [Atom(at_gro, at_itp)
                      for at_gro, at_itp in zip(mol._molecule_gro,  # type: ignore
                                                mol._molecule_itp)]
        # mol._init_atoms_bonds()
        return mol

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
        return self._atoms.index(atom)

    def hash2atom(self, number: int) -> 'Atom':
        """
        Returns the atom with the corresponding number as hash.

        This method is duplicated in MoleculeItp to avoid problems in getattr.

        Parameters
        ----------
        number : int
            The hash of the atom to access.

        Returns
        -------
        atom : AtomItp
            The atom with the number.

        """
        return self[self.hash2index(number)]

    @property
    def bonds_distance(self) -> Dict[int, List[Tuple[int, float]]]:
        """
        dict of int: list of tuple(int, float): A complex data structure
            that collect the information of the bond distances. The key of the
            property corresponds to the atom index in the molecule. The value
            is a list with tuples. For each tuple, the first value corresponds
            with the index of the bonded atom and the second is the length of
            the bond. This property is used in the alignment process in gaddle
            maps.
        """
        bond_info: Dict[int, List[Tuple[int, float]]] = defaultdict(list)
        for index, atom in enumerate(self):
            for hash_to in atom.bonds:
                index_to = self.hash2index(hash_to)
                atom_to = self[index_to]
                distance = euclidean(atom.position, atom_to.position)
                bond_info[index].append((index_to, distance))
        return dict(bond_info)

    @property
    def info(self) -> InfoDict:
        """
        dict : A dictionary with the needed information to restore the
            Molecule.
        """
        gro_info = [at.gro_line() for at in self.molecule_gro]  # type: ignore
        info = InfoDict()
        info.molecule_itp = os.path.abspath(self.molecule_itp.fitp)
        info.molecule_gro = gro_info
        return info

    @classmethod
    def from_info(cls, info_dict: InfoDict) -> 'Molecule':
        """
        Loads a molecule from the info porperty.

        Parameters
        ----------
        info_dict : A dictionary with the needed information to restore the
            object.

        Returns
        -------
        molecule: Molecule
            The loaded object.
        """
        
        mol_itp = MoleculeItp(info_dict.molecule_itp)
        if len(mol_itp.resnames) == 1:
            mol_gro = MoleculeGro([AtomGro(line)
                                   for line in info_dict.molecule_gro])
        else:
            atom_lines: List[AtomGro] = []
            molecules = []
            resname = ''
            for line in info_dict.molecule_gro:
                if line[1] != resname:
                    if atom_lines:
                        molecules.append(MoleculeGro(atom_lines))
                    atom_lines = []
                    resname = line[1]
                atom_lines.append(AtomGro(line))
            molecules.append(MoleculeGro(atom_lines))
            mol_gro = MacromoleculeGro(*molecules)

        return Molecule(mol_gro, mol_itp)

    @classmethod
    def from_gro_itp(cls, fgro: str, fitp: str) -> 'Molecule':
        """
        Loads the molecule from gro and itp file.

        Parameters
        ----------
        fgro : str
            The file name with the gro.
        fitp : str
            The file name with the itp.

        Returns
        -------
        molecule : Molecule
            The initialized molecule.

        """
        from . import System
        syst = System(fgro, fitp)
        if len(syst) == 1:
            return syst[0]
        raise IOError('The input gro file has more than one molecule.')


def _molecule_itp_and_residues_match(molecule_itp: MoleculeItp,
                                     residues: List[Residue]) -> bool:
    if len(molecule_itp) != sum(len(res) for res in residues):
        return False
    index = 0
    for res in residues:
        for atom in res:  # type: ignore
            if atom.residname != molecule_itp[index]:
                return False
            index += 1
    return True
