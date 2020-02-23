#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This module contains MoleculeGro and AtomGro objects which allows to load
atom and molecule information from .gro lines.

"""

import warnings
import collections
import re
import numpy as np
from typing import Tuple, Union, List, Optional

from ..parsers import GroFile
from . import GeneralAtom, GeneralMolecule


class MoleculeGro(GeneralMolecule):
    """
    Contains the information of a molecule atoms.

    Parameters
    ----------
    atoms : list of Atoms
        A list with the atoms of the molecule

    """

    def __init__(self, atoms):
        super(MoleculeGro, self).__init__()
        self._atoms_gro = atoms

    def __add__(self, other: Union['MacromoleculeGro',
                                   'AtomGro',
                                   'MoleculeGro']) -> Union['MacromoleculeGro',
                                                            'MoleculeGro']:
        if isinstance(other, MacromoleculeGro):
            # Only can be added if the last molecule of self has different
            # residname.
            if self.residname != other.residname[0]:
                molecules = [self] + other.molecules
                return MacromoleculeGro(*molecules)
            else:
                raise ValueError(('To sum a MoleculeGro with a MacromoleculeGro'
                                  ' the residname ({}) must be different to the'
                                  ' first residname of MacromoleculeGro ({}).'
                                  '').format(self.residname,
                                             other.residname[0]))
        elif isinstance(other, MoleculeGro):
            # Only can be added if they have the same residname
            if other.residname == self.residname:
                return MoleculeGro(self.atoms+other.atoms)
            return MacromoleculeGro(self, other)
        elif isinstance(other, AtomGro):
            # Only can be added if they have the same residname
            if other.parent_residname == self.residname:
                return MoleculeGro(self.atoms+[other.copy()])
            else:
                raise ValueError(('To sum an atom and a molecule '
                                  ' they has to have the same resid and '
                                  'resname. Given: {} to sum with'
                                  '{}').format(other.parent_residname,
                                               self.residname))
        else:
            raise TypeError(('It is not possible to sum a MoleculeGro and a {}'
                             '').format(type(other)))

    def __radd__(self, other):
        if not other:
            return self
        return self + other

    def __getitem__(self, index) -> 'AtomGro':
        return self._atoms_gro[index]

    def __len__(self) -> int:
        return len(self._atoms_gro)

    def __str__(self):
        string = 'Molecule of {} with resid {}.'.format(self.resname,
                                                        self.resid)
        return string

    __repr__ = __str__

    def __eq__(self, element):
        if isinstance(element, MoleculeGro):
            return super(MoleculeGro, self).__eq__(element)
        return False

    def __ne__(self, element):
        return not self == element

    @property
    def atoms(self) -> List['AtomGro']:
        """
        :obj:`list` of :obj:`AtomGro` : List with the atoms of the molecule.
        """
        return [atom.copy() for atom in self._atoms_gro]

    @property
    def atoms_positions(self) -> np.ndarray:
        """
        numpy.ndarray : An array with the atoms positions.
        """
        return np.array([a.position for a in self])

    @atoms_positions.setter
    def atoms_positions(self, new_positions: np.ndarray):
        for atom, pos in zip(self, new_positions):
            atom.position = pos

    @property
    def atoms_velocities(self) -> np.ndarray:
        """
        numpy.ndarray : An array with the atoms velocities.
        """
        return np.array([a.velocity for a in self])

    @atoms_velocities.setter
    def atoms_velocities(self, new_velocities: np.ndarray):
        if new_velocities is None:
            for atom in self:
                atom.velocity = None
            return 
        for atom, vel in zip(self, new_velocities):
            atom.velocity = vel

    @property
    def atoms_ids(self) -> List[int]:
        """
        list : A list with the ids of the atoms molecule.
        """
        return [at.atomid for at in self]

    @atoms_ids.setter
    def atoms_ids(self, new_ids: List[int]):
        if len(self) != len(new_ids):
            raise IndexError('The new ids must have the same length as self.')
        for atom, _id in zip(self, new_ids):
            atom.atomid = _id

    @property
    def geometric_center(self) -> np.ndarray:
        """
        :obj:`numpy.ndarray` (3) :
        Coordinates of the geometric center of the molecule.
        """
        return np.mean(self.atoms_positions, axis=0)

    @property
    def x(self) -> float:
        """
        The x coordinate of the geometric center
        """
        return self.geometric_center[0]

    @property
    def y(self) -> float:
        """
        The y coordinate of the geometric center
        """
        return self.geometric_center[1]

    @property
    def z(self) -> float:
        """
        The z coordinate of the geometric center
        """
        return self.geometric_center[2]

    @property
    def resname(self) -> str:
        """
        :obj:`string` : Resname of the molecule.
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

    @property
    def resid(self) -> Optional[int]:
        """
        :obj:`int` : Resid of the molecule.
        """
        try:
            return self[0].resid
        except IndexError:
            # In case of empty molecules
            return None

    @resid.setter
    def resid(self, value: int):
        for atom in self:
            atom.resid = value

    @property
    def residname(self) -> str:
        """
        string : An idtentier of the molecule (resid+name)
        """
        return '{}{}'.format(self.resid, self.resname)

    def rotate(self, rotation_matrix: np.ndarray):
        """
        Rotate the molecule around its center of mass with a given rotation
        matrix.

        Parameters
        ----------
        rotation_matrix : numpy.ndarray
            3x3 array with the rotation matrix.
            Advice: construct the matrix with the help of "rotation_matrix"
            function placed in compyna.

        """
        # Se desplaza el centro de masas de la molécula al origen se
        # rota y se desplaza
        com = self.geometric_center
        atoms_pos = self.atoms_positions - com
        new_pos = np.dot(atoms_pos, np.transpose(rotation_matrix)) + com
        self.atoms_positions = new_pos

    def remove_atom(self, atom: 'AtomGro'):
        """
        Removes a given atom from the molecule.

        Parameters
        ----------
        atom : AtomGro
             The atom we want to remove from the molecule.
        """
        self._atoms_gro.remove(atom)

    def move_to(self, new_position: np.ndarray):
        """
        Moves the molecule geometric_center to new_position.

        Parameters
        ----------
        new_position : numpy.ndarray(3)
            An array with the new position coordinates.

        """
        displacement = new_position - self.geometric_center
        self.move(displacement)

    def move(self, displacement: np.ndarray):
        """
        Displaces the molecule given a displacement vector.

        Parameters
        ----------
        displacement : numpy.ndarray(3)
            An array with the displacement vector.

        """
        self.atoms_positions = self.atoms_positions + displacement

    def plot(self, colors=None, atom_size=.1, fig=None, show=False):
        """
        This method makes a mayavi plot of the molecule.

        This method does not return nothing, remember to run mlab.show() to
        display the result.

        Parameters
        ----------
        colors : dictionary, optional
            One can specify the colors of the different atoms. Available colors
            are listed in function colors_codes(). If None a default
            set of colors will be used.
        atom_size : float
            The size for the sphere representing the atoms.
        fig : mayavi.figure (optional)
            The mayavi scene where the plot will be made.

        Returns
        -------
        fig : mayavi.figure
            The mayavi scene with the plot.

        Examples
        --------
        >>> colors = {'H':'white', 'C':'black', 'N':'blue'}

        """

        from .analyze.representations import (asign_atoms_color, color_codes,
                                              DEFAULT_COLORMAP)
        from mayavi import mlab

        fig = mlab.figure(fig)

        if colors is None:
            colors = DEFAULT_COLORMAP
        # Se completan los colores
        colors = asign_atoms_color([at.atomname for at in self], colors)
        # Se representan los átomos
        codes = color_codes()
        for atom in self:
            coord = atom.position
            mlab.points3d(
                *coord,
                scale_factor=atom_size,
                resolution=20,
                color=codes[colors[atom.atomname]],
                scale_mode='none')
        if show:
            mlab.show()
        return fig

    def copy(self) -> 'MoleculeGro':
        """
        Returns a copy of the molecule.

        Returns
        -------
        molecule : MoleculeGro
            The copy of the molecule.

        """
        mol = self.__new__(self.__class__)
        mol._atoms_gro = self.atoms
        return mol

    def write_gro(self, fout: str):
        """
        Writes a .gro file with the molecule conformation.

        Parameters
        ----------
        fout : str
            The path with the file to write the file.

        """
        with GroFile(fout, 'w') as fgro:
            for atom in self:
                fgro.writeline(atom.gro_line())

    def update_from_molecule_itp(self, mitp):
        """
        Modifies the atoms information to match the itp.

        Parameters
        ----------
        mtip : MoleculeItp
            The molecule to match

        Raises
        ------
        ValueError
            If number of atoms in self and in mitp does not match.

        """
        if len(mitp) != len(self):
            raise ValueError(('There gro molecule with resname {} in '
                              'are different from that in the itp.'
                              '').format(self.resname))
        for at_gro, at_itp in zip(self, mitp):
            at_gro.atomname = at_itp.atomname

    def distance_to(self, mol, box_vects=None, inv=False) -> float:
        """
        Returns the distance between self and mol.

        mol can be a molecules instance or a 3D vector.

        Parameters
        ----------
        mol : MoleculeGro or numpy.ndarray
            The molecule or a point to compute the distance.
        box_vects : numpy.ndarray
            The box vectors to apply pbcs.
        inv : bool
            If it is True, box_vects are considered as the inverse to afford
            calc.

        Returns
        -------
        distance : float
            The euclidean distance.

        """
        if isinstance(mol, MoleculeGro):
            mol = mol.geometric_center
        vect = mol - self.geometric_center

        if box_vects is not None:
            if not inv:
                box_vects = np.linalg.inv(box_vects)
            vect = vect.dot(box_vects)
            vect -= np.round(vect)
            vect = vect.dot(box_vects)

        return np.linalg.norm(vect)

    def align_to(self, atom_sel):
        """
        Aligns self to an atom selection.

        Calculates the self.atoms_positions coordinates in an axis base
        determined by the atoms with the same name as in atom_sel and then
        that positions are projected in the same base determined by atom_sel.

        Parameters
        ----------
        atom_sel : MDAnalysis.atomgroup
            An atom selection with 2, 3, 4 or 5 atoms (the number of positions
            allowed in calculate_ba

        """
        from .analyze.orientation import calcule_base
        sel_pos, sel_names = zip(*[(atom.position, atom.name)
                                   for atom in atom_sel])
        mol_pos = [atom.position for name in sel_names
                   for atom in self if atom.atomname == name]

        sel_base, sel_app = calcule_base(sel_pos)
        mol_base, mol_app = calcule_base(mol_pos)
        self.move(sel_app - mol_app)
        conversion_matrix = np.matmul(np.linalg.inv(sel_base), mol_base)
        self.atoms_positions = np.matmul(conversion_matrix,
                                         self.atoms_positions.T).T

    @property
    def distance_to_zero(self) -> float:
        """
        float : The distance between the geometric_center and (0, 0, 0)
            without applying pbcs.
        """
        return np.linalg.norm(self.geometric_center)


class MacromoleculeGro(MoleculeGro):
    """
    A class to wrap molecules that have more than one residue.

    This class is based based on MoleculeGro. The difference
    is that now the property resname returns a list with the found resnames.

    Parameters
    ----------
    *molecules_gro : *MoleculeGro
        A list with the MoleculeGro objects that conforms the macromolecule.

    Raises
    ------
    IOError
        If less than two molecule molecules are input or every molecules has
        the same resname.
    TypeError
        If one of the input molecule is not instance of MoleculeGro.

    """

    def __init__(self, *molecules_gro):

        # Input validation
        if len(molecules_gro) < 2:
            raise IOError(('To initialize a MacromoleculeGro more than two '
                           'MoleculeGro instance.'))

        if len(set(mol.residname for mol in molecules_gro)) == 1:
            raise IOError(('The input molecules have the same resname and resid'
                           ', try to generate a MoleculeGro instance.'))

        self._resnames_list = resnames_list = []
        self._molecules = []
        atoms = []
        for mol in molecules_gro:
            if not isinstance(mol, MoleculeGro):
                raise TypeError(('The input molecules should be instance of '
                                 'MoleculeGro not {}.').format(type(mol)))
            resnames_list.append(mol.resname)
            self._molecules.append(mol.copy())
            for atom in mol:
                atoms.append(atom)

        super(MacromoleculeGro, self).__init__(atoms)

    def __add__(self, other) -> 'MacromoleculeGro':
        if isinstance(other, MacromoleculeGro):
            if other.residname[-1] != self.residname[0]:
                molecules = self.molecules + other.molecules
                return MacromoleculeGro(*molecules)
            else:
                raise ValueError(('To sum two MacromoleculeGro, the last'
                                  ' residname ({}) of the first must be '
                                  'different to the first residname of the '
                                  'second ({}).').format(other.residname,
                                                         self.residname[-1]))
        elif isinstance(other, MoleculeGro):
            # Only can be added if the last molecule of self has different
            # residname.
            if other.residname != self.residname[-1]:
                molecules = self.molecules + [other]
                return MacromoleculeGro(*molecules)
            else:
                raise ValueError(('To sum a MacromoleculeGro with a MoleculeGro'
                                  ' the residname ({}) must be different to the'
                                  ' last residname of MacromoleculeGro ({}).'
                                  '').format(other.residname,
                                             self.residname[-1]))
        else:
            raise TypeError(('It is not possible to sum a MacromoleculeGro and'
                             ' a {}').format(type(other)))

    @property
    def resname(self) -> List[str]:
        """
        :obj:`list` : Resname list of the molecule.
        """
        return self._resnames_list

    @property
    def resnames(self) -> List[str]:
        """
        :obj:`list` : Resname list of the molecule.
        """
        return self.resname

    @resname.setter
    def resname(self, new_resname):
        if isinstance(new_resname, (list, tuple)):
            if len(new_resname) != self._resnames_list:
                raise IOError(('The number of input resnames must be {}'
                               '').format(len(new_resname)))
            for mol, res in zip(self._molecules, new_resname):
                res = GroFile.validate_string(res)
                mol.resname = res
            self._resnames_list = new_resname[:]
        elif isinstance(new_resname, str):
            res = GroFile.validate_string(new_resname)
            for mol in self._molecules:
                mol.resnmae = res
        else:
            raise TypeError('Wrong resname value for MacromoleculeGro.')

    @property
    def resid(self) -> List[int]:
        """
        :obj:`int` : Resid of the molecule.
        """
        return [mol.resid for mol in self._molecules]

    @resid.setter
    def resid(self, value):
        if isinstance(value, collections.abc.Iterable):
            for mol, resid in zip(self._molecules, value):
                mol.resid = resid
        elif isinstance(value, int):
            for mol in self._molecules:
                mol.resid = value
        else:
            raise TypeError('Wrong resid value for MacromoleculeGro.')

    @property
    def residname(self) -> List[str]:
        """
        string : An idtentier of the molecule (resid+name)
        """
        return ['{}{}'.format(resid, res)
                for resid, res in zip(self.resid, self.resname)]

    @property
    def molecules(self) -> List['MoleculeGro']:
        """
        list of MoleculeGro : A list with the simple molecules conforming self.
        """
        return [mol.copy() for mol in self._molecules]

    def copy(self) -> 'MacromoleculeGro':
        """
        Returns a copy of the macromolecule.

        Returns
        -------
        molecule : MacromoleculeGro
            The copy of the molecule.

        """
        mol = self.__new__(self.__class__)
        mol._molecules = self.molecules
        mol._resnames_list = self._resnames_list
        mol._atoms_gro = self.atoms
        return mol


class AtomGro(GeneralAtom):
    """
    Returns an object with the information of a .gro line corresponding to an
    atom.

    Parameters
    ----------
    parsed_gro_line : list
        A list returned by the read_gro function when a correct formatted .gro
        line is input.

    """
    # pylint: disable=too-few-public-methods

    def __init__(self, parsed_gro_line):
        super(AtomGro, self).__init__()
        (self.resid,
         self.resname,
         self.atomname,
         self.atomid,
         cordx,
         cordy,
         cordz) = parsed_gro_line[:7]

        self.position = np.array([cordx, cordy, cordz])
        velocities = parsed_gro_line[7:10]
        if not velocities:
            self.velocity = None
        else:
            self.velocity = np.array(velocities)

    def __add__(self, other) -> 'MoleculeGro':
        if isinstance(other, AtomGro):
            if other.parent_residname == self.parent_residname:
                return MoleculeGro([self, other])
            else:
                raise ValueError(('To sum an atom to another atom'
                                  ' they has to have the same resid and '
                                  'resname. Given: {} to sum with'
                                  '{}').format(other.parent_residname,
                                               self.parent_residname))
        else:
            # Try the commutative add to see if the other object
            # has the addition defined
            return other + self

    def __str__(self):
        string = 'Atom {} of molecule {} with number {}.'.format(self.atomname,
                                                                 self.resname,
                                                                 self.resid)
        return string

    __repr__ = __str__

    def gro_line(self, parsed=True) -> List: 
        """
        Returns the gro line corresponding to the atom.

        Parameters
        ----------
        parsed : bool
            If True, the line is returned as read_gro function output.
            Else, formated line is returned.

        """
        elements = [
            self.resid,
            self.resname,
            self.atomname,
            self.atomid,
        ]
        elements += list(self.position)
        if self.velocity is not None:
            elements += list(self.velocity)
        if parsed:
            return elements
        return GroFile.parse_atomlist(elements)

    def copy(self) -> 'AtomGro':
        """
        Returns a copy of the atom.

        Returns
        -------
        new_atom : AtomGro
            The copied atom.

        """
        atom = self.__new__(self.__class__)
        atom.resid = self.resid
        atom.resname = self.resname
        atom.atomname = self.atomname
        atom.atomid = self.atomid
        atom.position = self.position.copy()
        if self.velocity is None:
            atom.velocity = None
        else:
            atom.velocity = self.velocity.copy()
        return atom

    @property
    def parent_residname(self) -> str:
        """
        string : The corresponding molecule hash.
        """
        return '{}{}'.format(self.resid, self.resname)

    @property
    def element(self) -> str:
        """
        string : The type of the atom (element).
        """
        element = re.findall(r'([A-Za-z]+)', self.atomname)
        if element:
            return element[0]
        else:
            raise IOError(('Wrong format for atomname ({}). The element can'
                           ' not be parsed.'.format(self.atomname)))
