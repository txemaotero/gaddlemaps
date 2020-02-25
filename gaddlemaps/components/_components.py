# -*- coding: utf-8 -*-
'''
This module contains Atom and Molecule objects which connect information from
both gro and itp files.
'''

import warnings
import os
from collections import defaultdict
from typing import Optional, Any, List, Dict, Tuple, Union
from scipy.spatial.distance import euclidean
from . import (AtomGro, AtomItp, GeneralAtom, GeneralMolecule, MoleculeItp,
               MoleculeGro, MacromoleculeGro)


class Atom(GeneralAtom):
    """
    An atom class that wraps the AtomItp and AtomGro classes.

    Properties and methods are the same as both AtomItp and AtomGro.

    Parameters
    ----------
    atom_gro : AtomGro
        The AtomGro object.
    atom_itp : AtomItp
        The AtomItp object.

    Raises
    ------
    IOError
        If the atom_gro and atom_itp do not correspond to the same atom.
    TypeError
        If the inputs are not instance of the corresponding class.

    """

    def __init__(self, atom_gro: AtomGro, atom_itp: AtomItp):
        super(Atom, self).__init__()
        self._atom_gro: AtomGro = None
        self._atom_itp: AtomItp = None
        self.atom_gro = atom_gro
        self.atom_itp = atom_itp

    @property
    def atom_gro(self) -> AtomGro:
        """
        AtomGro : The input AtomGro object
        """
        return self._atom_gro

    @atom_gro.setter
    def atom_gro(self, new_atom: AtomGro):
        if not self._atom_itp is None:
            if not GeneralAtom.are_the_same_atom(new_atom, self._atom_itp):
                msg = 'The input atoms do not correspond to the same atom.'
                raise IOError(msg)
        if not isinstance(new_atom, AtomGro):
            raise TypeError('atom_gro input have to be an AtomGro instance.')
        self._atom_gro = new_atom

    @property
    def atom_itp(self) -> AtomItp:
        """
        AtomItp : The input AtomItp object
        """
        return self._atom_itp

    @atom_itp.setter
    def atom_itp(self, new_atom: AtomItp):
        if self._atom_gro is not None:
            if not GeneralAtom.are_the_same_atom(new_atom, self._atom_gro):
                msg = 'The input atoms do not correspond to the same atom.'
                raise IOError(msg)
        if not isinstance(new_atom, AtomItp):
            raise TypeError('atom_itp input have to be an AtomItp instance.')
        self._atom_itp = new_atom

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._atom_gro, attr):
            return getattr(self._atom_gro, attr)
        elif hasattr(self._atom_itp, attr):
            return getattr(self._atom_itp, attr)
        else:
            raise AttributeError('Atom object has no attribute {}'.format(attr))

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_atom_itp', '_atom_gro']:
            super(Atom, self).__setattr__(attr, value)
        # Special cases that are present in both atoms
        elif attr == 'resname':
            setattr(self._atom_itp, attr, value)
            setattr(self._atom_gro, attr, value)
        elif attr in super(Atom, self).__dir__():
            super(Atom, self).__setattr__(attr, value)
        elif hasattr(self._atom_itp, attr):
            setattr(self._atom_itp, attr, value)
        elif hasattr(self._atom_gro, attr):
            setattr(self._atom_gro, attr, value)
        else:
            super(Atom, self).__setattr__(attr, value)

    def __str__(self) -> str:
        string = 'Atom {} of molecule {} with number {}.'.format(self.atomname,
                                                                 self.resname,
                                                                 self.resid)
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
        at = self.__new__(self.__class__)
        at._atom_gro = self._atom_gro.copy()
        at._atom_itp = self._atom_itp
        return at

    def __eq__(self, atom: Any) -> bool:
        cond = super(Atom, self).__eq__(atom)
        if cond:
            if isinstance(atom, (Atom, AtomItp)):
                return self.number == atom.number
            return cond
        return False


class Molecule(GeneralMolecule):
    """
    Loads a molecule itp and the .gro information.

    This class wraps all the features of both MoleculeItp and MoleculeGro
    objects. Saves a copy of the molecule_gro and the original molecule_itp.

    Parameters
    ----------
    molecule_itp : MoleculeItp
        The object with the .itp information of the molecule.
    molecule_gro : MoleculeGro
        The object with the .gro information of the molecule.

    Raises
    ------
    TypeError
        If the input are instances of wrong type.
    IOError
        If molecule_gro and molecule_itp represent different molecules.

    """
    def __init__(self, molecule_gro: MoleculeGro, molecule_itp: MoleculeItp):
        self._molecule_itp: MoleculeItp = None
        self._molecule_gro: MoleculeGro = None
        self.molecule_itp = molecule_itp
        self.molecule_gro = molecule_gro
        self._atoms = [Atom(at_gro, at_itp)
                       for at_gro, at_itp in zip(self._molecule_gro,
                                                 self._molecule_itp)]

    @property
    def molecule_itp(self) -> MoleculeItp:
        """
        MoleculeItp : The object with the .itp information of the molecule.
        """
        return self._molecule_itp

    @molecule_itp.setter
    def molecule_itp(self, new_molecule: MoleculeItp):
        if not self._molecule_gro is None:
            if not self.are_the_same_molecule(self.molecule_gro, new_molecule):
                raise IOError(('The input molecule_itp does not match with the'
                               'loaded molecule_gro.'))
        if not isinstance(new_molecule, MoleculeItp):
            msg = ('The first input should be an instance of MoleculeItp not'
                   ' {}').format(type(new_molecule))
            raise TypeError(msg)
        self._molecule_itp = new_molecule

    @property
    def molecule_gro(self) ->  MoleculeGro:
        """
        MoleculeGro : The object with the .gro information of the molecule.
        """
        return self._molecule_gro

    @molecule_gro.setter
    def molecule_gro(self, new_molecule: MoleculeGro):
        if not self._molecule_itp is None:
            if not self.are_the_same_molecule(self.molecule_itp, new_molecule):
                raise IOError(('The input molecule_gro does not match with the'
                               'loaded molecule_itp.'))
        if not isinstance(new_molecule, MoleculeGro):
            msg = ('molecule_gro should be an instance of MoleculeGro not'
                   ' {}').format(type(new_molecule))
            raise TypeError(msg)
        self._molecule_gro = new_molecule.copy()

    def __getitem__(self, index: int) -> 'Molecule':
        return self._atoms[index]

    def __len__(self) -> int:
        return len(self._atoms)

    def __iter__(self) -> 'Atom':
        for atom in self._atoms:
            yield atom

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._molecule_itp, attr):
            return getattr(self._molecule_itp, attr)
        elif hasattr(self._molecule_gro, attr):
            return getattr(self._molecule_gro, attr)
        else:
            raise AttributeError(('Molecule object has no attribute {}'
                                  '').format(attr))

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_molecule_itp', '_molecule_gro']:
            super(Molecule, self).__setattr__(attr, value)
        elif attr in super(Molecule, self).__dir__():
            super(Molecule, self).__setattr__(attr, value)
        elif hasattr(self._molecule_itp, attr):
            setattr(self._molecule_itp, attr, value)
        elif hasattr(self._molecule_gro, attr):
            setattr(self._molecule_gro, attr, value)
        else:
            super(Molecule, self).__setattr__(attr, value)

    def __str__(self) -> str:
        string = 'Molecule of {}.'.format(self.name)
        return string

    __repr__ = __str__

    def __eq__(self, molecule: Any) -> bool:
        if isinstance(molecule, Molecule):
            return super(Molecule, self).__eq__(molecule)
        return False

    def __ne__(self, molecule: Any) -> bool:
        return not self == molecule

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Molecule, self).__dir__())
        _dir.update(dir(self._molecule_gro))
        _dir.update(dir(self._molecule_itp))
        return list(_dir)

    @property
    def atoms(self) -> List['Atom']:
        """
        list of Atom: List with the copies of the atoms of the molecule.
        """
        return [atom.copy() for atom in self._atoms]

    @property
    def resname(self) -> str:
        """
        string: Resname of the molecule.
        """
        return self[0].resname

    @resname.setter
    def resname(self, new_resname: str):
        if isinstance(new_resname, str):
            if len(new_resname) > 5:
                old_resname = new_resname
                new_resname = new_resname[:5]
                warn_text = ('The input resname has more than 5 character ({})'
                             ', this is modified to be {}.').format(old_resname,
                                                                    new_resname)
                warnings.warn(warn_text, UserWarning)
            for atom in self:
                atom.resname = new_resname
        else:
            raise TypeError('Resname must be a string.')

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
            mol._molecule_gro = self._molecule_gro.copy()
        else:
            if GeneralMolecule.are_the_same_molecule(self._molecule_gro,
                                                     new_molecule_gro):
                mol._molecule_gro = new_molecule_gro.copy()
            else:
                raise IOError(('The input molecule_gro can not replace the old'
                               ' one, please try to create a new instance.'))
        mol._atoms = [Atom(at_gro, at_itp)
                      for at_gro, at_itp in zip(mol._molecule_gro,
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
        bond_info = defaultdict(list)
        for index, atom in enumerate(self):
            for hash_to in atom.bonds:
                index_to = self.hash2index(hash_to)
                atom_to = self[index_to]
                distance = euclidean(atom.position, atom_to.position)
                bond_info[index].append((index_to, distance))
        return dict(bond_info)

    @property
    def info(self) -> Dict[str, Union[str, List[Union[str, int, float]]]]:
        """
        dict : A dictionary with the needed information to restore the
            Molecule.
        """
        gro_info = [at.gro_line() for at in self.molecule_gro]
        info = {
            'molecule_itp': os.path.abspath(self.molecule_itp.fitp),
            'molecule_gro': gro_info,
        }
        return info

    @classmethod
    def from_info(cls, info_dict: Dict[str, Union[str, List[Union[str, int, float]]]]) -> 'Molecule':
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
        mol_itp = MoleculeItp(info_dict['molecule_itp'])
        if len(mol_itp.resnames) == 1:
            mol_gro = MoleculeGro([AtomGro(line)
                                   for line in info_dict['molecule_gro']])
        else:
            atom_lines = []
            molecules = []
            resname = ''
            for line in info_dict['molecule_gro']:
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
