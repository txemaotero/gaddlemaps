# -*- coding: utf-8 -*-
"""
This module contains MoleculeItp and AtomItp objects which allows to load
atom and molecule information from .itp files.
"""

import warnings
from itertools import groupby

from typing import List, Dict, Any, Union, Tuple, Set

from ..parsers import ItpLineAtom, ItpFile, ItpLineBonds
from . import GeneralAtom, GeneralMolecule


class MoleculeItp:
    """
    Loads molecules from .itp topology file.

    This class behaves like a list of atoms which has bonds defined.

    Parameters
    ----------
    fitp : string
        The path  to the .itp file with the molecule information.

    Raises
    ------
    IOError
        If the input itp file miss information.

    """
    def __init__(self, fitp: str):
        self.fitp = fitp
        itp_file = ItpFile(fitp)

        if 'moleculetype' not in itp_file:
            raise IOError('The input itp must have "moleculetype" section.')
        if 'atoms' not in itp_file:
            raise IOError('The input itp must have "atoms" section.')

        self.name = ''
        self._atoms_itp: List['AtomItp'] = []
        self._number_to_index: Dict[int, int] = {}

        self._init_name(itp_file)
        self._init_atoms(itp_file)
        self._init_atoms_bonds(itp_file)

    def __getitem__(self, index: int) -> 'AtomItp':
        return self._atoms_itp[index]

    def __len__(self) -> int:
        return len(self._atoms_itp)

    def __eq__(self, element: Any) -> bool:
        if isinstance(element, MoleculeItp):
            if element.name == self.name and len(self) == len(element):
                return all(at1 == at2 for at1, at2 in zip(self, element)) # type: ignore
        return False

    def __ne__(self, element: Any) -> bool:
        return not self == element

    def __str__(self) -> str:
        string = f'Molecule of {self.name} from itp.'
        return string

    __repr__ = __str__

    def _init_name(self, itp_file: ItpFile):
        for line in itp_file['moleculetype']:
            if line.content:
                self.name = line.name
        if not self.name:
            raise IOError('There is not molecule name in "moleculetype" sec.')

    def _init_atoms(self, itp_file: ItpFile):
        atoms_sec = itp_file['atoms']
        index = 0
        for atom_line in atoms_sec:
            if atom_line.content:
                atom = AtomItp(atom_line)
                self._atoms_itp.append(atom)
                self._number_to_index[hash(atom)] = index
                index += 1
        if not self._atoms_itp:
            raise IOError('There are not atoms in the atoms section.')

    def _init_atoms_bonds(self, itp_file: ItpFile):
        # If there is more than one atom, there has to be bonds.
        atoms = self._atoms_itp
        condition_constraints = 'constraints' in itp_file
        condition_bonds = 'bonds' in itp_file
        if len(atoms) > 1:
            if not (condition_constraints or condition_bonds):
                # Special case of TP5 (only settles and not bonds)
                if not (('settles' in itp_file) or
                        ('virtual_sites2' in itp_file) or
                        ('virtual_sites3' in itp_file) or
                        ('virtual_sites4' in itp_file)):
                    raise IOError('The input file has to provide atoms bonds')
            bonds: List[ItpLineBonds] = []
            if condition_constraints:
                bonds += list(itp_file['constraints'])
            if condition_bonds:
                bonds += list(itp_file['bonds'])
            index_parser = {atom.number: atom for atom in atoms}
            for bond in bonds:
                index_parser[bond.atom_from].connect(index_parser[bond.atom_to])
        else:
            if condition_bonds or condition_constraints:
                warnings.warn('{}: There are bonds specifications but the '
                              'molecule has only one atom.'.format(self.name))

    @property
    def resnames(self) -> List[str]:
        """
        list of str: Resnames of the atoms without consecutive repetitions.
        """
        # TODO: Ver si se puede usar este property para que devuelva una lista
        # con todos los nombres de residues (aunque sea repetido) y luego que
        # al seter le puedas pasar un str (todos con el mismo) o una lista con
        # todos los residuos aunque sea repetidos. Pensar cual sería la mejor
        # forma.
        tot_resnames = ('{:5}{}'.format(atom.resname, atom.resnr)
                        for atom in self)  # type: ignore
        return [x[0][:5].strip() for x in groupby(tot_resnames)]

    @property
    def resname_len_list(self) -> List[Tuple[str, int]]:
        """
        list of tuple(str, int) : 
            [(resname_1, number_of_atoms_with_resname_1),
            (resname_2, number_of_atoms_with_resname_2), ...]
        """
        tot_resnames = ('{:5}{}'.format(atom.resname, atom.resnr)
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

    def hash2atom(self, number: int) -> 'AtomItp':
        """
        Returns the atom with the corresponding number as hash.

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

    def hash2index(self, atom_hash: int) -> int:
        """
        Returns the index of an atom by its hash.

        Parameters
        ----------
        atom_hash : int
            The hash of the atom (it usually starts in 1).

        Returns
        -------
        atom_index : int
            The index of the atom in the molecule (it usually starts in 0).
        """
        return self._number_to_index[atom_hash]

    def index2hash(self, atom_index: int) -> int:
        """
        The reverse operation done by hash2index.
        """
        return hash(self[atom_index])

    def index(self, atom: 'AtomItp') -> int:
        """
        Returns the index of the atom in the molecule.

        Parameters
        ----------
        atom : AtomItp
            The atom to find the index.

        Returns
        -------
        index : int
            The index of the atom in the molecule.

        """
        return self._atoms_itp.index(atom)


class AtomItp:
    """
    Contains the information of a .itp line corresponding to an atom.

    The instances of this class have access to all the information accessible
    through the ItpLineAtom class.

    Parameters
    ----------
    itp_line_atom : ItpLineAtom or string
        A parsed atom line from the itp (see parse module) or a string
        corresponding to a line in the itp with atomic information.

    Attributes
    ----------
    bonds: set of int
        A set with the hash of the atoms that are connected to self.

    """
    def __init__(self, itp_line_atom: Union[str, ItpLineAtom]):
        if isinstance(itp_line_atom, ItpLineAtom):
            self._itp_line_atom = itp_line_atom
        elif isinstance(itp_line_atom, str):
            self._itp_line_atom = ItpLineAtom(itp_line_atom)
        self.bonds: Set[int] = set()
        # Overwrite mutable attributes
        self.resname = self._itp_line_atom.resname
        self.atomname = self._itp_line_atom.atomname

    def __getattr__(self, attr: str) -> Any:
        return getattr(self._itp_line_atom, attr)

    def __repr__(self) -> str:
        return f'Itp atom of {self.atomname} of molecule {self.resname}'

    __str__ = __repr__

    def __hash__(self) -> int:
        """
        Number in the itp line
        """
        return self.number

    def __eq__(self, atom: 'AtomItp') -> bool:  # type: ignore
        if isinstance(atom, AtomItp):
            condition = (
                (self.number == atom.number) and
                (self.resname == atom.resname) and
                (self.atomname == atom.atomname)
            ) 
            return condition
        return False

    @property
    def residname(self) -> str:
        """
        string: An identifier of the residue (resid+name)
        """
        return '{}{}'.format(self.resid, self.resname)

    def connect(self, atom: 'AtomItp'):
        """
        Connects self with other atom setting the bond.

        Parameters
        ----------
        atom : AtomItp
            Atom to connect.

        """
        self.bonds.add(hash(atom))
        atom.bonds.add(hash(self))

    def closest_atoms(self, natoms: int = 2) -> List[int]:
        """
        Returns a list natoms bonded atoms to self.

        If more than natoms are bonded self, the natoms  with lower id_num are
        returned.

        Parameters
        ----------
        natoms : integer
            The number of atoms to return.

        Returns
        -------
        bonded_atoms : list of AtomItp
            The list with the natoms bonded atoms.
        """
        return sorted(self.bonds)[:natoms]
