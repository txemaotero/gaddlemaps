# -*- coding: utf-8 -*-
"""
This module defines the System and SystemGro classes.
"""

import os
from itertools import groupby, product, islice
from collections import Counter
from typing import List, Tuple, Union, Generator, Set, Mapping, Dict

import numpy as np

from . import Residue, AtomGro, MoleculeItp, Molecule
from ..parsers import GroFile


class System(object):
    """
    Class to manage simulation systems.

    A System object is formed by Molecule objects. Only the molecules
    corresponding to the input fitps will be loaded.

    Parameters
    ----------
    fgro : string
        Gromacs file with the system information.
    *fitps : string
        Paths with the itp files to load molecules.

    """
    def __init__(self, fgro: str, *fitps: str):
        self.system_gro = SystemGro(fgro)
        self.different_molecules: List[Molecule] = []
        # (ind_mol, ind_gro_start, ammount)
        self._molecules_ordered: List[Tuple[int, int, int]] = []
        self._available_mgro_ordered = np.array(list(self.system_gro.molecules_info_ordered_all))
        for fitp in fitps:
            self.add_fitp(fitp)

    def __str__(self) -> str:
        string = 'Simulation system with:\n\n'
        string += '\n'.join(["{:6}: {}".format(k, v)
                             for k, v in sorted(self.composition.items())])
        return string

    __repr__ = __str__

    def __iter__(self) -> Molecule:
        for index, gro_start, gro_end in self._molecules_ordered_all_gen():
            m_gro = sum(self.system_gro[gro_start:gro_end])
            mol = self.different_molecules[index].copy(m_gro)
            yield mol

    def __getitem__(self, index: int) -> Molecule:
        if isinstance(index, slice):
            molecules = []
            info = list(self._molecules_ordered_all_gen())
            for info in islice(self._molecules_ordered_all_gen(), index.start,
                               index.stop, index.step):
                itp_index, gro_start, gro_end = info
                m_gro = sum(self.system_gro[gro_start:gro_end])
                mol = self.different_molecules[itp_index].copy(m_gro)
                molecules.append(mol)
            return molecules
        for i, info in enumerate(self._molecules_ordered_all_gen()):
            if i == index:
                itp_index, gro_start, gro_end = info
                m_gro = sum(self.system_gro[gro_start:gro_end])
                mol = self.different_molecules[itp_index].copy(m_gro)
                return mol
        raise IndexError('Molecule index out of range')

    def __len__(self) -> int:
        return sum(elem[2] for elem in self._molecules_ordered)

    def _check_index_in_available_mgro(self, index_array: np.ndarray, 
                                       mol_itp: MoleculeItp) -> int:
        start_index = None
        len_array = len(index_array)
        for mgro_index, mgro_pk in enumerate(self._available_mgro_ordered):
            if mgro_pk == index_array[0]:
                if (self._available_mgro_ordered[mgro_index:mgro_index+len_array] == index_array).all():
                    start_index = mgro_index
                    break
        if start_index is None:
            raise IOError(('The sequence of residues in found in '
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

    def _molecules_ordered_all_gen(self) -> Generator:
        for index, gro_start, ammount in self._molecules_ordered:
            len_mol = len(self.different_molecules[index].resnames)
            for i in range(ammount):
                yield (index, gro_start+i*len_mol, gro_start+(i+1)*len_mol)

    def add_fitp(self, fitp: str):
        """
        Adds and identifies the molecule from the fitp to the system.
        """
        mol_itp = MoleculeItp(fitp)
        self.add_molecule_itp(mol_itp)

    def add_molecule_itp(self, mol_itp: MoleculeItp):
        """
        Adds a molecule to the system and find it in the gro file.
        """
        # Check if all the residues in the itp are in the gro file
        index_mol_gro = []
        gro_mols_resnames = self.system_gro.molecules_resname_len_index
        for resname_len_itp in mol_itp.resname_len_list:
            if resname_len_itp not in gro_mols_resnames:
                raise IOError(f'The molecule {mol_itp.name} is not in the gro file')
            index_mol_gro.append(gro_mols_resnames[resname_len_itp])
        # Check if the resnames appear in the same order in the gro file.
        index_mol_gro = np.array(index_mol_gro)
        start_index = self._check_index_in_available_mgro(index_mol_gro,
                                                          mol_itp)
        # Try to init the molecule
        mol_gro = sum(self.system_gro[start_index:start_index+len(index_mol_gro)])
        molecule = Molecule(mol_gro, mol_itp)
        mol_index = len(self.different_molecules)
        self.different_molecules.append(molecule)
        self._find_all_molecules_and_replace(index_mol_gro, mol_index,
                                             start_index)
        # Sort the molecules in ordered_molecules
        self._molecules_ordered.sort(key=lambda x: x[1])

    @property
    def composition(self) -> Mapping[str, int]:
        """
        Counter of str: int : For each resname (key), how many molecules
            there are (value).
        """
        composition = Counter()
        for index, _, ammount in self._molecules_ordered:
            composition[self.different_molecules[index].name] += ammount
        return composition

    @property
    def info(self) -> Dict[str, Union[str, List[str]]]:
        """
        dict : A dictionary with the needed information to restore the
            system.

        Example:

        info = {
            'system_gro': '/home/system.gro',
            'mols_itps': ['/home/mol1.itp', '/home/mol2.itp'],
        }

        """
        info = {
            'system_gro': self.fgro,
            'mols_itps': self.fitps,
        }
        return info

    @property
    def fgro(self) -> str:
        return self.system_gro.fgro

    @property
    def fitps(self) -> List[str]:
        """
        list of string : A list with the paths to the itp files of the loaded
            molecules.
        """
        return [os.path.abspath(mol.fitp) for mol in self.different_molecules]

    @classmethod
    def from_info(cls, info_dict: Dict[str, Union[str, List[str]]]) -> 'System':
        """
        Builds the system from the information returned by info property.

        Parameters
        ----------
        info_dict : A dictionary with the needed information to restore the
            system.

            Example:

            info = {
                'system_gro': '/home/system.gro',
                'mols_itps': ['/home/mol1.itp', '/home/mol2.itp'],
            }

        Returns
        -------
        system : System
            The loaded system.
        """
        return System(info_dict['system_gro'], *info_dict['mols_itps'])


class SystemGro(object):
    """
    Class to work with the information in gro files.

    Basically this class acts as a list of Residue objects. Only one Residue
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
        self._open_fgro = GroFile(fgro)
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

    def __getitem__(self, index: int) -> Union[Residue, List[Residue]]:
        if isinstance(index, slice):
            molecules = []
            for _, start, len_mol in islice(self._molecules_ordered_all_gen(),
                                            index.start, index.stop,
                                            index.step):
                self._open_fgro.seek_atom(start)
                molecules.append(Residue([AtomGro(next(self._open_fgro))
                                          for _ in range(len_mol)]))
            return molecules
        for i, info in enumerate(self._molecules_ordered_all_gen()):
            if i == index:
                _, start, len_mol = info
                self._open_fgro.seek_atom(start)
                return Residue([AtomGro(next(self._open_fgro))
                                for _ in range(len_mol)])
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

    @box_matrix.setter
    def box_matrix(self, new_matrix: np.ndarray):
        self._open_fgro.box_matrix = new_matrix

    @property
    def comment_line(self) -> str:
        """
        str: The comment line in the gro file.
        """
        return self._open_fgro.comment

    @comment_line.setter
    def comment_line(self, new_comment: str):
        self._open_fgro.comment = new_comment

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
