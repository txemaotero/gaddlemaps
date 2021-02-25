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
This submodule defines the System and SystemGro classes.
"""

import typing
from collections import Counter
from more_itertools import islice_extended, last
from typing import Dict, Generator, List, Mapping, Tuple, overload

import numpy as np

from ..parsers import open_coordinate_file
from . import AtomGro, Molecule, MoleculeTop, Residue


class System:
    """
    Class to manage simulation systems.

    A System object is formed by Molecule objects. Only the molecules
    corresponding to the input ftops will be loaded.

    Parameters
    ----------
    fgro : string
        Gromacs file with the system information.
    *ftops : string
        Paths with the files with the bonds information to load molecules.

    Raises
    ------
    IOError
        If one of the topology files do not match with any molecule in the
        system.

    """
    def __init__(self, fgro: str, *ftops: str):
        self.system_gro = SystemGro(fgro)
        self.different_molecules: List[Molecule] = []
        # (ind_mol, ind_gro_start, ammount)
        self._molecules_ordered: List[List[int]] = []
        self._available_mgro_ordered = np.array(list(self.system_gro.molecules_info_ordered_all))
        for ftop in ftops:
            self.add_ftop(ftop)

    def __str__(self) -> str:
        if not self.composition:
            return 'Simulation system with no loaded molecules.'
        string = 'Simulation system with:\n\n'
        string += '\n'.join(["{:6}: {}".format(k, v)
                             for k, v in sorted(self.composition.items())])
        return string

    __repr__ = __str__

    def __iter__(self) -> Generator[Molecule, None, None]:
        for index, gro_start, gro_end in self._molecules_ordered_all_gen():
            residues = self.system_gro[gro_start:gro_end]  # type: ignore
            mol = self.different_molecules[index].copy(residues)
            yield mol

    @overload
    def __getitem__(self, index: int) -> Molecule:
        ...
    @overload
    def __getitem__(self, index: slice) -> List[Molecule]:
        ...
    def __getitem__(self, index):
        try:
            if isinstance(index, slice):
                molecules = []
                info = list(self._molecules_ordered_all_gen())
                for info in islice_extended(self._molecules_ordered_all_gen(),
                                            index.start, index.stop,
                                            index.step):
                    itp_index, gro_start, gro_end = info
                    residues = self.system_gro[gro_start:gro_end]
                    mol = self.different_molecules[itp_index].copy(residues)
                    molecules.append(mol)
                return molecules
            elif isinstance(index, int):
                if index == -1:
                    info = last(self._molecules_ordered_all_gen())
                else:
                    info = next(islice_extended(self._molecules_ordered_all_gen(),
                                                index, index+1))
                itp_index, gro_start, gro_end = info
                residues = self.system_gro[gro_start:gro_end]
                mol = self.different_molecules[itp_index].copy(residues)
                return mol
            else:
                raise TypeError(('System indices must be integers or slices,'
                                 f' not {type(index)}.'))
        except StopIteration:
            raise IndexError('Molecule index out of range')

    def __len__(self) -> int:
        return sum(elem[2] for elem in self._molecules_ordered)

    def _check_index_in_available_mgro(self, index_array: np.ndarray,
                                       mol_itp: MoleculeTop) -> int:
        start_index = None
        len_array = len(index_array)
        for mgro_index, mgro_pk in enumerate(self._available_mgro_ordered):
            if mgro_pk == index_array[0]:
                if (self._available_mgro_ordered[mgro_index:mgro_index+len_array] == index_array).all():
                    start_index = mgro_index
                    break
        if start_index is None:
            raise IOError(('The sequence of residues found for '
                           f'{mol_itp.name} is not found in the gro file, so '
                           'the molecule is not recognized in the system.'))
        return start_index

    def _find_all_molecules_and_replace(self, index_mol_gro: np.ndarray,
                                        mol_index: int,
                                        start_index: int = 0):
        av_gro = self._available_mgro_ordered
        l_index_mol = len(index_mol_gro)
        new_block = True
        while (start_index+l_index_mol) <= len(av_gro):
            if (av_gro[start_index:start_index+l_index_mol] == index_mol_gro).all():
                if new_block:
                    self._molecules_ordered.append([mol_index, start_index, 1])
                    new_block = False
                else:
                    self._molecules_ordered[-1][2] += 1
                av_gro[start_index:start_index+l_index_mol] = -1
                start_index += l_index_mol
            else:
                new_block = True
                start_index += 1

    def _molecules_ordered_all_gen(self) -> Generator[Tuple[int, int, int],
                                                      None, None]:
        for index, gro_start, ammount in self._molecules_ordered:
            len_mol = len(self.different_molecules[index].resnames)
            for i in range(ammount):
                yield (index, gro_start+i*len_mol, gro_start+(i+1)*len_mol)

    def add_ftop(self, ftop: str):
        """
        Adds and identifies the molecule from the ftop to the system.
        """
        self.add_molecule_top(MoleculeTop(ftop))

    def add_molecule_top(self, mol_top: MoleculeTop):
        """
        Adds a molecule to the system and find it in the gro file.
        """
        # Check if all the residues in the itp are in the gro file
        index_mol_gro = []
        gro_mols_resnames = self.system_gro.molecules_resname_len_index
        for resname_len_top in mol_top.resname_len_list:
            if resname_len_top not in gro_mols_resnames:
                raise IOError(f'The molecule {mol_top.name} is not in the gro file')
            index_mol_gro.append(gro_mols_resnames[resname_len_top])
        # Check if the resnames appear in the same order in the gro file.
        index_mol_gro = np.array(index_mol_gro)
        start_index = self._check_index_in_available_mgro(index_mol_gro,
                                                          mol_top)
        # Try to init the molecule
        residues = self.system_gro[start_index:start_index + len(index_mol_gro)]
        molecule = Molecule(mol_top, residues)
        mol_index = len(self.different_molecules)
        self.different_molecules.append(molecule)
        self._find_all_molecules_and_replace(index_mol_gro, mol_index,
                                             start_index)
        # Sort the molecules in ordered_molecules
        self._molecules_ordered.sort(key=lambda x: x[1])

    @property
    def composition(self) -> typing.Counter[str]:
        """
        Counter of str: int : For each molecule name (key), how many
            molecules there are (value).
        """
        composition: typing.Counter[str] = Counter()
        for index, _, ammount in self._molecules_ordered:
            composition[self.different_molecules[index].name] += ammount
        return composition

    @property
    def fgro(self) -> str:
        return self.system_gro.fgro


class SystemGro:
    """
    Class to work with the information in gro files.

    Basically this class acts like a list of Residue objects. Only one Residue
    instance of each type is loaded to afford memory for large systems. The
    positions of the rest of the residues are storage and they are generated
    when requested.

    Parameters
    ----------
    fgro : string
        Gromacs file with the system information.

    """
    def __init__(self, fgro: str):
        self.fgro = fgro
        self._open_fgro = open_coordinate_file(fgro)
        self.different_molecules: List[Residue] = []
        self._molecules_pk: Dict[Tuple[str, int], int] = {}   # Dictionary with {(resname, len(mol)): index}
        self._molecules_ordered: List[int] = []
        self._parse_gro()

    def __str__(self) -> str:
        string = 'Simulation system with:\n\n'
        string += '\n'.join(["{:6}: {}".format(k, v)
                             for k, v in sorted(self.composition.items())])
        return string

    __repr__ = __str__

    def __iter__(self) -> Generator[Residue, None, None]:
        for _, start, len_mol in self._molecules_ordered_all_gen():
            self._open_fgro.seek_atom(start)
            yield Residue([AtomGro(next(self._open_fgro))
                           for _ in range(len_mol)])

    def __del__(self):
        self._open_fgro.close()

    @overload
    def __getitem__(self, index: int) -> Residue:
        ...
    @overload
    def __getitem__(self, index: slice) -> List[Residue]:
        ...
    def __getitem__(self, index):
        try:
            if isinstance(index, slice):
                molecules = []
                for _, start, len_mol in islice_extended(self._molecules_ordered_all_gen(),
                                                        index.start, index.stop,
                                                        index.step):
                    self._open_fgro.seek_atom(start)
                    molecules.append(Residue([AtomGro(next(self._open_fgro))
                                              for _ in range(len_mol)]))
                return molecules
            if isinstance(index, int):
                if index == -1:
                    info = last(self._molecules_ordered_all_gen())
                else:
                    info = next(islice_extended(self._molecules_ordered_all_gen(),
                                                index, index+1))
                _, start, len_mol = info
                self._open_fgro.seek_atom(start)
                return Residue([AtomGro(next(self._open_fgro))
                                for _ in range(len_mol)])
            else:
                raise TypeError(('SystemGro indices must be integers or slices,'
                                 f' not {type(index)}.'))
        except StopIteration:
            raise IndexError('Residue index out of range')

    def __len__(self) -> int:
        return sum(elem[1] for elem in self._pk_ammount_ordered_gen())

    def _parse_gro(self):
        current_residue = [AtomGro(next(self._open_fgro))]
        prev_atom_residname = current_residue[0].residname
        for line in self._open_fgro:
            atom = AtomGro(line)
            if atom.residname == prev_atom_residname:
                current_residue.append(atom)
            else:
                self._add_residue_init(Residue(current_residue))
                current_residue = [atom]
                prev_atom_residname = current_residue[0].residname
        self._add_residue_init(Residue(current_residue))

    def _add_residue_init(self, residue: Residue):
        key = (residue.resname, len(residue))
        if residue not in self.different_molecules:
            self.different_molecules.append(residue)
            index = len(self.different_molecules) - 1
            self._molecules_pk[key] = index
        index = self._molecules_pk[key]

        if not self._molecules_ordered or self._molecules_ordered[-2] != index:
            self._molecules_ordered += [index, 1]
        else:
            self._molecules_ordered[-1] += 1

    def _pk_ammount_ordered_gen(self) -> Generator[Tuple[int, int], None, None]:
        for i in range(0, len(self._molecules_ordered), 2):
            yield (self._molecules_ordered[i], self._molecules_ordered[i+1])

    def _molecules_ordered_all_gen(self) -> Generator[Tuple[int, int, int], None, None]:
        start_atom = 0
        for index, ammount in self._pk_ammount_ordered_gen():
            len_mol = len(self.different_molecules[index])
            for _ in range(ammount):
                yield (index, start_atom, len_mol)
                start_atom += len_mol

    @property
    def n_atoms(self) -> int:
        """
        int : The number of atoms in the system.
        """
        return self._open_fgro.natoms

    @property
    def box_matrix(self) -> np.ndarray:
        """
        numpy.ndarray (3,3) : the 3 lattice vectors of the box.
        """
        return self._open_fgro.box_matrix

    @property
    def comment_line(self) -> str:
        """
        str: The comment line in the gro file.
        """
        return self._open_fgro.comment

    @property
    def molecules_info_ordered_all(self) -> Generator[int, None, None]:
        """
        generator of int: Returns the index of the molecules in the
            different_molecules attribute in order of appearance.
        """
        for index, ammount in self._pk_ammount_ordered_gen():
            for _ in range(ammount):
                yield index

    @property
    def molecules_resname_len_index(self) -> Dict[Tuple[str, int], int]:
        """
        dict of tuple (string, int) to int: The keys of the dictionary are
            tuples with the residue name and the number of atoms of the
            different residues in the system and the values are its index in
            the list of different molecules.
        """
        return self._molecules_pk

    @property
    def composition(self) -> Mapping[str, int]:
        """
        Counter of str: int : For each resname (key), how many molecules
            there are (value).
        """
        composition: Mapping[str, int] = Counter()
        for index, ammount in self._pk_ammount_ordered_gen():
            composition[self.different_molecules[index].resname] += ammount  # type: ignore
        return composition
