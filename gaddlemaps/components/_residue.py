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
This submodule defines a Residue class that will constitute molecules.
"""

import re
import warnings
from typing import (TYPE_CHECKING, Any, Generator, List, Optional, Union,
                    overload)

import numpy as np

from ..parsers import open_coordinate_file, GroLine, GroFile

if TYPE_CHECKING:
    from . import MoleculeTop


class Residue:
    """
    A class with the information of a residue in a .gro file.

    This class creates objects that are enumerations of AtomGro instances. It
    has methods to manipulate atoms positions maintaining the shape of the
    residue. This class have to be initialized with non empty list of atoms.

    You can add two residues if the atoms from both residues have the same
    residue name and number. This can also be done with just one atom and the
    same check will be done.

    Parameters
    ----------
    atoms : list of AtomGro
        A list with the atoms of the residue.

    Raises
    ------
    ValueError : If the the input list of atoms is empty or they have different
        residue number or name.
    """

    def __init__(self, atoms: List['AtomGro']):
        residnames = set(atom.residname for atom in atoms)
        if not residnames:
            raise ValueError(('Empty list of atoms are not allowed to '
                              'Residue initialization'))
        if len(residnames) > 1:
            raise ValueError(('The atoms of a Residue must have the same'
                              ' residue name and number. Found residues: '
                              f'{residnames}.'))
        self._atoms_gro = atoms


    @overload
    def __getitem__(self, index: int) -> 'AtomGro':
        ...
    @overload
    def __getitem__(self, index: slice) -> List['AtomGro']:
        ...
    def __getitem__(self, index):
        return self._atoms_gro[index]

    def __iter__(self) -> Generator['AtomGro', None, None]:
        for atom in self._atoms_gro:
            yield atom

    def __len__(self) -> int:
        return len(self._atoms_gro)

    def __str__(self) -> str:
        return f'Residue {self.resname} with number {self.resid}.'

    __repr__ = __str__

    def __eq__(self, element: Any) -> bool:
        """
        Two residues are equal if all their atoms are equal.
        """
        if isinstance(element, Residue) and (len(self) == len(element)):
            return all(at1 == at2 for at1, at2 in zip(element, self))
        return False

    def __ne__(self, element: Any) -> bool:
        return not self == element

    def __add__(self, other: Union['Residue', 'AtomGro']) -> 'Residue':
        if isinstance(other, Residue):
            if self.residname != other.residname:
                raise ValueError(('To add two Residues they must have the same'
                                  ' residue name and number. Given residues:'
                                  f' {self.residname}, {other.residname}.'))
            return Residue(self.atoms + other.atoms)
        elif isinstance(other, AtomGro):
            # Only can be added if they have the same residname
            if other.residname == self.residname:
                return Residue(self.atoms + [other.copy()])
            else:
                raise ValueError(('To add an AtomGro and a Residue '
                                  ' they must have the same residue name and'
                                  ' number. Given atom residue: '
                                  f'{other.residname}. Current residue: '
                                  f'{self.residname}.'))
        else:
            raise TypeError(('It is not possible to add a Residue and a '
                             f'{type(other)}'))

    def __radd__(self, other: Union['Residue', 'AtomGro']) -> 'Residue':
        if not other:
            return self.copy()
        return self + other

    @property
    def atoms(self) -> List['AtomGro']:
        """
        list of AtomGro: List with a copy of the atoms of the residue.
        """
        return [atom.copy() for atom in self._atoms_gro]

    @property
    def atoms_positions(self) -> np.ndarray:
        """
        numpy.ndarray((N, 3)) : An array with the atoms positions.
        """
        return np.array([a.position for a in self])

    @atoms_positions.setter
    def atoms_positions(self, new_positions: np.ndarray):
        if new_positions.shape != (len(self), 3):
            raise ValueError(('The new positions must be an array of shape '
                              f'({len(self)}, 3)'))
        for atom, pos in zip(self, new_positions):
            atom.position = pos

    @property
    def atoms_velocities(self) -> Optional[np.ndarray]:
        """
        numpy.ndarray((N, 3)) or None : An array with the atoms velocities.
            If one of the atoms has no velocity this returns None.
        """
        vels = []
        for atom in self:
            if atom.velocity is None:
                return None
            vels.append(atom.velocity)
        return np.array(vels)

    @atoms_velocities.setter
    def atoms_velocities(self, new_velocities: Optional[np.ndarray]):
        if new_velocities is None:
            for atom in self:
                atom.velocity = None
            return
        if new_velocities.shape != (len(self), 3):
            raise ValueError(('The new velocities must be an array of shape '
                              f'({len(self)}, 3)'))
        for atom, vel in zip(self, new_velocities):
            atom.velocity = vel

    @property
    def atoms_ids(self) -> List[int]:
        """
        list of int: A list with the ids of the atoms of the residue.
        """
        return [at.atomid for at in self]

    @atoms_ids.setter
    def atoms_ids(self, new_ids: List[int]):
        if len(self) != len(new_ids):
            raise IndexError('The new ids must have the same length as self.')
        if not all(isinstance(i, int) for i in new_ids):
            raise TypeError(f'atomids must be integers, given: {new_ids}')
        for atom, _id in zip(self, new_ids):
            atom.atomid = _id

    @property
    def geometric_center(self) -> np.ndarray:
        """
        numpy.ndarray(3): Coordinates of the geometric center of the residue.
        """
        return np.mean(self.atoms_positions, axis=0)

    @property
    def x(self) -> float:
        """
        float: The x coordinate of the geometric center of the residue.
        """
        return self.geometric_center[0]

    @property
    def y(self) -> float:
        """
        float: The y coordinate of the geometric center of the residue.
        """
        return self.geometric_center[1]

    @property
    def z(self) -> float:
        """
        float: The z coordinate of the geometric center of the residue.
        """
        return self.geometric_center[2]

    @property
    def resname(self) -> str:
        """
        string: Resname of the residue.
        """
        return self[0].resname

    @resname.setter
    def resname(self, new_resname: str):
        if isinstance(new_resname, str):
            if len(new_resname) > 5:
                old_resname = new_resname
                new_resname = new_resname[:5]
                warn_text = ('The input resname has more than 5 character '
                             f' ({old_resname}). This has been modified to'
                             f' {new_resname}.')
                warnings.warn(warn_text, UserWarning)
            for atom in self:  # type: ignore
                atom.resname = new_resname
        else:
            raise TypeError('Resname must be a string.')

    @property
    def resid(self) -> int:
        """
        int: Residue number of the residue.
        """
        return self[0].resid

    @resid.setter
    def resid(self, value: int):
        for atom in self:
            atom.resid = value

    @property
    def residname(self) -> str:
        """
        string: An identifier of the residue (resid+name)
        """
        return '{}{}'.format(self.resid, self.resname)

    def rotate(self, rotation_matrix: np.ndarray):
        """
        Rotate the residue around its center of mass with a given rotation
        matrix.

        Parameters
        ----------
        rotation_matrix : numpy.ndarray
            3x3 array with the rotation matrix.
            ADVICE: create the matrix with the help of "rotation_matrix"
            function placed in the root of the package.

        """
        com = self.geometric_center
        atoms_pos = self.atoms_positions - com
        new_pos = np.dot(atoms_pos, np.transpose(rotation_matrix)) + com
        self.atoms_positions = new_pos

    def remove_atom(self, atom: 'AtomGro'):
        """
        Removes a given atom from the residue.

        Parameters
        ----------
        atom : AtomGro
             The atom you want to remove from the residue.
        """
        self._atoms_gro.remove(atom)

    def move_to(self, new_position: np.ndarray):
        """
        Moves the residue geometric_center to new_position.

        Parameters
        ----------
        new_position : numpy.ndarray(3)
            An array with the new position coordinates.

        """
        displacement = new_position - self.geometric_center
        self.move(displacement)

    def move(self, displacement: np.ndarray):
        """
        Moves the residue a given displacement vector.

        Parameters
        ----------
        displacement : numpy.ndarray(3)
            An array with the displacement vector.

        """
        self.atoms_positions = self.atoms_positions + displacement

    def copy(self) -> 'Residue':
        """
        Returns a copy of the residue.

        Returns
        -------
        residue : Residue
            The copy of the residue.

        """
        return Residue(self.atoms)

    def write_gro(self, fout: str):
        """
        Writes a .gro file with the residue conformation.

        Parameters
        ----------
        fout : str
            The path with the file to write the information.

        """
        with open_coordinate_file(fout, 'w') as fgro:
            for atom in self:
                fgro.writeline(atom.gro_line())

    def update_from_molecule_top(self, mtop: "MoleculeTop"):
        """
        Modifies the Residue atoms name to match the mtop.

        This method is very useful when you have a miss-match between the
        atom names in the topology and gro files. This will modify the Residue
        atoms names to match the names in the topology. Make sure that all the
        atoms in the topology are in the Residue.

        Parameters
        ----------
        mtop : MoleculeTop
            The molecule to match.

        Raises
        ------
        ValueError
            If number of atoms in self and in mtop does not match.

        """
        if len(mtop) != len(self):
            raise ValueError((f'The Residue with name {self.resname} '
                              f'missmatch the itp molecule {mtop.name}.'))
        for at_gro, at_itp in zip(self, mtop):
            at_gro.name = at_itp.name

    def distance_to(self, residue: Union['Residue', np.ndarray],
                    box_vects: np.ndarray = None,
                    inv: bool = False) -> float:
        """
        Returns the distance between self and residue.

        residue can be a Residue instance or a 3D vector.

        Parameters
        ----------
        residue : Residue or numpy.ndarray
            The residue or a point to compute the distance.
        box_vects : numpy.ndarray
            The box vectors to apply periodic boundary conditions.
        inv : bool
            If it is True, box_vects are considered as the inverse matrix of
            the actual box_vects for a better performance.

        Returns
        -------
        distance : float
            The euclidean distance.

        """
        if isinstance(residue, Residue):
            residue = residue.geometric_center
        vect = residue - self.geometric_center

        if box_vects is not None:
            if not inv:
                box_vects = np.linalg.inv(box_vects)
            vect = vect.dot(box_vects)
            vect -= np.round(vect)
            vect = vect.dot(box_vects)

        return np.linalg.norm(vect)

    @property
    def distance_to_zero(self) -> float:
        """
        float : The distance between the geometric_center and (0, 0, 0)
            without applying periodic boundary conditions.
        """
        return np.linalg.norm(self.geometric_center)


class AtomGro:
    """
    It contains the information of a .gro line corresponding to an atom.

    You can add atoms with the same residue name and number to form a Residue
    instance.

    Parameters
    ----------
    parsed_gro_line : list
        A list returned by the GroFile.read_line method when a correct
        formatted .gro line is input.

    Attributes
    ----------
    resid : int
        Residue number of the atom.
    resname : str
        Residue name of the atom.
    name : str
        Name of the atom.
    atomid : str
        Atom index in the .gro file.
    position : np.ndarray(3)
        The coordinates of the atom.
    velocity : np.ndarray(3) or None
        The velocity of the atom if it is specified in the .gro line. If not
        this attribute is set to None.
    """

    def __init__(self, parsed_gro_line: GroLine):
        super(AtomGro, self).__init__()
        (self.resid,
         self.resname,
         self.name,
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

    def __add__(self, other: 'AtomGro') -> 'Residue':
        if isinstance(other, AtomGro):
            if other.residname == self.residname:
                return Residue([self, other])
            else:
                raise ValueError(('To sum an atom to another atom'
                                  ' they has to have the same resid and '
                                  'resname. Given: {} to sum with'
                                  '{}').format(other.residname,
                                               self.residname))
        else:
            # Try the commutative add to see if the other object
            # has the addition defined
            return other + self

    def __str__(self) -> str:
        string = 'Atom {} with residue {} and number {}.'.format(self.name,
                                                                 self.resname,
                                                                 self.resid)
        return string

    __repr__ = __str__

    def __eq__(self, element: Any) -> bool:
        """
        Two AtomGro are equal if they have the same name and the same residue
        name.
        """
        if isinstance(element, AtomGro):
            return ((self.resname == element.resname) and
                    (self.name == element.name))
        return False

    def __ne__(self, element: Any) -> bool:
        return not self == element

    def gro_line(self, parsed: bool = True) -> Union[str, List[Union[str, int, float]]]:
        """
        Returns the gro line corresponding to the atom.

        Parameters
        ----------
        parsed : bool, optional
            If True (default), the line is returned as FileGro.read_line
            method output. Else, line with the correct .gro format is
            returned.

        Returns
        -------
        gro_line: List of (str, int or float) or str
            The corresponding .gro line.

        """
        elements = [
            self.resid,
            self.resname,
            self.name,
            self.atomid,
        ]
        elements += list(self.position)
        if self.velocity is not None:
            elements += list(self.velocity)
        if parsed:
            return elements
        return GroFile.parse_atomlist(elements)  # type: ignore

    def copy(self) -> 'AtomGro':
        """
        Returns a copy of the atom.

        You can safely change the attributes of the returned atom without
        changing the original one.

        Returns
        -------
        new_atom : AtomGro
            The copied atom.

        """
        input_list = [self.resid, self.resname, self.name, self.atomid]
        input_list += list(self.position)
        if self.velocity is not None:
            input_list += list(self.velocity)
        return AtomGro(input_list)  # type: ignore

    @property
    def residname(self) -> str:
        """
        string : An identifier that contains the residue name and number. This
            should be unique for ear Residue object in a simulation system.
        """
        return '{}{}'.format(self.resid, self.resname)

    @property
    def element(self) -> str:
        """
        string : The element of the atom. It is obtained removing the non
            alphabetic character in the atom name.
        """
        element = re.findall(r'([A-Za-z]+)', self.name)
        if element:
            return element[0]
        else:
            raise IOError(('Wrong format for name ({}). The element can'
                           ' not be parsed.'.format(self.name)))
