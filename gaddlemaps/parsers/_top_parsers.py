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
This module defines functions to parse topology files and return a list of the
atoms already bonded.
"""

from typing import TYPE_CHECKING, List, Tuple, Dict, Optional, Type
import abc
import os.path

from . import ItpFile, ItpLineAtom, ItpLineBonds

if TYPE_CHECKING:
    from ..components import AtomTop


class TopologyParserManager:
    parsers: Dict[str, Type['TopologyParser']] = {}

    @classmethod
    def register(cls, parser: Type['TopologyParser']):
        if parser.EXTENSIONS:
            for extension in parser.EXTENSIONS:
                cls.parsers[extension] = parser

class TopologyParserRegistered(abc.ABCMeta):
    def __init__(self, name, bases, attrs):
        super().__init__(name, bases, attrs)
        TopologyParserManager.register(self)

    def __new__(metaclass, name, bases, attrs):
        return super().__new__(metaclass, name, bases, attrs)

class TopologyParser(metaclass=TopologyParserRegistered):

    EXTENSIONS: Optional[Tuple[str, ...]] = None

    @abc.abstractmethod
    def __init__(self, fit: str):
        super().__init__()

    @property
    @abc.abstractmethod
    def molecule_name(self) -> str:
        """
        Returns the name of the molecule in the file
        """
        return ""

    @property
    @abc.abstractmethod
    def atoms_info(self) -> List[Tuple[str, str, int]]:
        """
        Returns a list of tuples with atomic info. Each tuple must contain, the
        name of the residue in which the atom is included, the name of the atom,
        and an index of the atom. The index must be unique inside the molecule.
        """

        example = [
            ("residue_name", "atom_name", 1),
            ("residue_name", "atom_name_2", 2),
            ("residue_name_2", "atom_name_3", 3),
            ]
        return example


    @property
    @abc.abstractmethod
    def atoms_bonds(self) -> List[Tuple[int, int]]:
        """
        Returns a list of tuples. The is 1 tuple per bond in the molecule. Each
        tuple has 2 interger that correspond to the atomic indexes (as output by
        atoms_info) of the atoms that take part in the bond.
        """
        example = [
            (0, 1),
            (1, 2),
        ]
        return example

    @property
    def all_info(self) -> Tuple[str, List[Tuple[str, str, int]],
                                List[Tuple[int, int]]]:
        """
        Returns all the info in the following order:
        - molecule_name
        -atoms_info
        -atoms_bonds
        """
        return self.molecule_name, self.atoms_info, self.atoms_bonds


class ItpParser(TopologyParser):
    """
    Reads the itp file and returns the information to load a molecule.

    This function extracts the name of the molecule and the atoms and residues
    names. It also extract a list with the bonds.

    Parameters
    ----------
    fitp : str
        The itp file name.

    """


    EXTENSIONS = ("itp", "ITP")

    def __init__(self, fitp: str):
        itp_file = ItpFile(fitp)

        if 'moleculetype' not in itp_file:
            raise IOError('The input itp must have "moleculetype" section.')
        if 'atoms' not in itp_file:
            raise IOError('The input itp must have "atoms" section.')

        self._name = _itp_top_name(itp_file)
        self._atoms, self._bonds = _itp_top_atoms(itp_file)

    @property
    def molecule_name(self) -> str:
        """
        Returns the name of the molecule in the file
        """
        return self._name

    @property
    def atoms_info(self) -> List[Tuple[str, str, int]]:
        """
        Returns a list of tuples with atomic info. Each tuple must contain, the
        name of the residue in which the atom is included, the name of the atom,
        and an index of the atom. The index must be unique inside the molecule.
        """

        return self._atoms


    @property
    def atoms_bonds(self) -> List[Tuple[int, int]]:
        """
        Returns a list of tuples. The is 1 tuple per bond in the molecule. Each
        tuple has 2 interger that correspond to the atomic indexes (as output by
        atoms_info) of the atoms that take part in the bond.
        """

        return self._bonds


def read_topology(ftop: str,
                  file_format: Optional[str]=None) -> Tuple[str,
                                                            List[Tuple[str, str, int]],
                                                            List[Tuple[int, int]]]:
    """
    Reads a known topology file and returns the information to load a molecule.

    This function extracts the name of the molecule and the atoms and residues
    names. It also extract a list with the bonds.

    Parameters
    ----------
    ftop : str
        The topology file name.

    file_format: Optional[str]
        Force the file to be read as if it had the extension 'file_format'

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
    name = os.path.basename(ftop)
    if file_format is None:
        extension = name.split(".")[-1]
    else:
        extension = file_format

    if extension in TopologyParserManager.parsers:
        return TopologyParserManager.parsers[extension](ftop).all_info
    else:
        raise ValueError(f"No parser available for extension {extension}")


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
    for key in ('constraints', 'bonds', 'pairs'):
        for bond in itp_file.get(key, []):
            bonds.append((bond.atom_from, bond.atom_to))
    return bonds
