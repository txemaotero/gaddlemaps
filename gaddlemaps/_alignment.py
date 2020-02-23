#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
This module contains features to manage the alignment engine.
'''

import numpy as np

from .components import Molecule
from .parsers import GroFile
from . import ExchangeMap, minimize_molecules


class Alignment(object):
    """
    Class to manage the molecule alignment.

    Parameters
    ----------
    start : Molecule
        The molecule to align in the initial resolution.
    end : Molecule
        The molecule to align in the end resolution.

    Attributes
    ----------
    exchange_map : ExchangeMap
        The map to change from initial resolution to the end one. It has to be
        initialized with the init_exchange_map method.
    SIGMA_SCALE : float
        A number that modulates the molecule displacements.
    STEPS_FACTOR : int
        The factor to apply in the calculation of the number of steps in the
        Monte-Carlo alignment.

    """

    SIGMA_SCALE = 0.5
    STEPS_FACTOR = 5000

    def __init__(self, start=None, end=None):
        super(Alignment, self).__init__()
        self._start = None
        self._end = None
        self.start = start
        self.end = end
        self.exchange_map = None

    @property
    def start(self):
        """
        Molecule : The molecule in the initial resolution.

        Setter Raises
        -------------
        TypeError
            If the molecule is not instance of Molecule
        ValueError
            If the molecule does not match with the end.
        """
        return self._start

    @start.setter
    def start(self, molecule):
        if molecule is None:
            self._start = None
            return
        if not isinstance(molecule, Molecule):
            raise TypeError(('The start attribute has to be an instance of'
                             ' Molecule, not {}.'.format(type(molecule))))
        if (self._end is None) or (self._start is None):
            self._start = molecule.copy()
        else:
            if molecule == self._start:
                self._start = molecule.copy()
            else:
                raise ValueError(('The start molecule do not match with the'
                                  ' end one. Pleas check that both molecules'
                                  ' are the same.'))

    @property
    def end(self):
        """
        Molecule : The molecule in the initial resolution.

        Setter Raises
        -------------
        TypeError
            If the molecule is not instance of Molecule
        ValueError
            If the molecule does not match with the start.
        """
        return self._end

    @end.setter
    def end(self, molecule):
        if molecule is None:
            self._end = None
            return
        if not isinstance(molecule, Molecule):
            raise TypeError(('The end attribute has to be an instance of'
                             ' Molecule, not {}.'.format(type(molecule))))
        if (self._start is None) or (self._end is None):
            self._end = molecule.copy()
        else:
            if molecule == self._end:
                self._end = molecule.copy()
            else:
                raise ValueError(('The end molecule do not match with the'
                                  ' end one. Pleas check that both molecules'
                                  ' are the same.'))

    def align_molecules(self, restrictions=None, deformation_types=None,
                        ignore_hydrogens=True):
        """
        Starts the alignment engine to find the optimal overlap between molecs.

        NOTE: If None is input as restrictions and the molecules to align have
        multiple residues, they will be guessed by residue matching (see
        guess_protein_restrains function).

        Parameters
        ----------
        restrictions : list of tuple of int, optional
            A list of tuples with pairs of atom numbers corresponding to start
            and end atoms molecules. The align will be performed privileging
            configurations where those atoms are close. By default is set to [].

            Example:
            >>> restrictions = [(1, 3), (4, 5)]

            IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
            (STARTS IN 0).

        deformation_types : tuple of int, optional
            Specifies the type of the minimization. Possible options:
                0 : Translation
                1 : Rotation
                2 : Individual atom move
            If it is None, all the possibilities are chosen.
        ignore_hydrogens : bool, optional
            If True, hydrogen atoms will not be included in the minimization
            of the distances. This will be only applied to the molecule which
            is not moved ind the alignment engine.

        """

        if restrictions is None:
            # Check for multiple residues
            restrictions = []
            if len(self.start.resnames) > 1:
                try:
                    restrictions = guess_protein_restrains(self.start,
                                                           self.end)
                except IOError:
                    pass

        if deformation_types is None:
            if len(self.start) == 1 or len(self.end) == 1:
                deformation_types = (0, )
            else:
                deformation_types = (0, 1, 2)

        # Move molecules to share geometric center
        self.start.move_to(self.end.geometric_center)
        # Identify the smallest molecule
        if len(self.start) < len(self.end):
            molecules = [self.end, self.start]
            restrictions = [i[::-1] for i in restrictions]
        else:
            molecules = [self.start, self.end]
        molecules = sorted([self.start, self.end], key=len, reverse=True)
        # Init the inputs to the backend
        if ignore_hydrogens:
            mol1_positions, restrictions = remove_hydrogens(molecules[0],
                                                            restrictions)
        else:
            mol1_positions = molecules[0].atoms_positions
        mol2_positions = molecules[1].atoms_positions
        mol2_bonds_info = molecules[1].bonds_distance
        mol2_com = molecules[1].geometric_center
        n_steps = self.STEPS_FACTOR*len(molecules[1])
        mol1_bonds_info = molecules[0].bonds_distance.values()
        translation_width = 2*min(t[1] for l in mol1_bonds_info for t in l)
        mol2_positions = minimize_molecules(mol1_positions, mol2_positions,
                                            mol2_com, self.SIGMA_SCALE,
                                            n_steps, restrictions,
                                            mol2_bonds_info, translation_width,
                                            deformation_types)
        # Change mol2 positions. Ensure that the changes are applied.
        if len(self.start) < len(self.end):
            self.start.atoms_positions = mol2_positions
        else:
            self.end.atoms_positions = mol2_positions

    def init_exchange_map(self, scale_factor=0.5):
        """
        Initializes the exchange map with the current molecules configuration.

        Parameters
        ----------
        scale_factor : float
            The compression factor to apply to mapped molecules.

        """
        self.exchange_map = ExchangeMap(self.start, self.end, scale_factor)

    def write_comparative_gro(self, fname=None):
        """
        Writes a .gro file with start and end molecules for check overlap.

        The molecules have START and END resnames respectively to ease the
        representation.

        Parameters
        ----------
        fname : string, optional
            The output .gro file name. If it is None, fname is set as:
            {name}_compare.gro, where name is the name of the molecule.

        """
        if fname is None:
            fname = '{}_compare.gro'.format(self.start.name)
        start = self.start.molecule_gro.copy()
        end = self.end.molecule_gro.copy()
        start.resname = 'START'
        end.resname = 'END'
        start.resid = 1
        end.resid = 2
        start.atoms_velocities = None
        end.atoms_velocities = None
        with GroFile(fname, 'w') as fgro:
            for atom in start:
                fgro.writeline(atom.gro_line())
            for atom in end:
                fgro.writeline(atom.gro_line())

    @property
    def info(self):
        """
        dict : A dictionary with the needed information to restore the
            Alignment object.
        """
        ex_map = (None if self.exchange_map is None
                  else self.exchange_map.scale_factor)
        mol_dic = lambda mol: mol if mol is None else mol.info
        info = {
            'start': mol_dic(self.start),
            'end': mol_dic(self.end),
            'exchange_map': ex_map,
        }
        return info

    @classmethod
    def from_info(cls, info_dict):
        """
        Builds the Alignment from the information returned by info property.

        Parameters
        ----------
        info_dict : A dictionary with the needed information to restore the
            object.

        Returns
        -------
        align : Alignment
            The loaded object.
        """
        ali = Alignment(start=Molecule.from_info(info_dict['start']),
                        end=Molecule.from_info(info_dict['end']))
        if info_dict['exchange_map'] is not None:
            ali.init_exchange_map(info_dict['exchange_map'])
        return ali


def remove_hydrogens(molecule, restrictions):
    """
    Returns positions of atoms that are not hydrogens and fixes restrictions.

    Parameters
    ----------
    molecule : Molecule
        The molecule to remove hydrogens.
    restrictions : list of tuple of int
        A list of tuples with pairs of atom numbers corresponding to start
        and end atoms molecules. The align will be performed privileging
        configurations where those atoms are close. By default is set to [].

        Example:
        >>> restrictions = [(1, 3), (4, 5)]

        IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE

    Returns
    -------
    new_positions : numpy.ndarray
        The positions of the atoms that are not hydrogens.
    new_restrictions : list of tuple of int
        The restrictions with the correct indexes.

    """
    positions = []
    index_1map = {}
    for index, atom in enumerate(molecule):
        if atom.element != 'H':
            positions.append(atom.position)
            index_1map[index] = len(positions) - 1
    new_restrictions = []
    for index_1, index_2 in restrictions:
        if index_1 in index_1map:
            new_restrictions.append((index_1map[index_1], index_2))
    return np.array(positions), new_restrictions


def guess_protein_restrains(mol1, mol2):
    """
    Guess restriains for molecules with multiple residues.

    Checks if mol1 and mol2 have the same resnames and create restrains pairing
    atoms with the same residues. This function only works with molecules with
    multiple residues.

    Parameters
    ----------
    mol1 : Molecule or MoleculeGro
        The first molecule to find the restrains.
    mol2 : Molecule or MoleculeGro
        The second molecule to find the restrains.

    Returns
    -------
    restrains : list of tuple of int
        A list of tuples with pairs of atom numbers corresponding to start
        and end atoms molecules. The align will be performed privileging
        configurations where those atoms are close. By default is set to [].

        Example:
        >>> restrictions = [(1, 3), (4, 5)]

        IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
        (STARTS IN 0).

    Raises
    ------
    IOError
        If the input molecules has different resnames.
    """
    if len(mol1.resnames) != len(mol2.resnames):
        raise IOError(('Input molecules must have the same number of residues:'
                       ' {}, {}.'.format(len(mol1.resnames),
                                         len(mol2.resnames))))
    if mol1.resnames != mol2.resnames:
        # Special case where resnames are similar
        for res1, res2 in zip(mol1.resnames, mol2.resnames):
            if (res1 in res2) or (res2 in res1):
                continue
            raise IOError(('Input molecules must have the same residues in the'
                           ' same order.'))
    restrains = []
    offset1 = 0
    offset2 = 0
    for mol_res1, mol_res2 in zip(mol1.molecules, mol2.molecules):
        restrains += guess_molecule_restrains(mol_res1, mol_res2, offset1,
                                              offset2)
        offset1 += len(mol_res1)
        offset2 += len(mol_res2)
    return restrains


def guess_molecule_restrains(mol1, mol2, offset1=0, offset2=0):
    """
    Guess restriains for molecules with one residue.

    Parameters
    ----------
    mol1 : Molecule or MoleculeGro
        The first molecule to find the restrains.
    mol2 : Molecule or MoleculeGro
        The second molecule to find the restrains.
    offset1 : int
        An offset to add to the atom index of mol1.
    offset2 : int
        An offset to add to the atom index of mol2.

    Returns
    -------
    restrains : list of tuple of int
        A list of tuples with pairs of atom numbers corresponding to start
        and end atoms molecules. The align will be performed privileging
        configurations where those atoms are close. By default is set to [].

        Example:
        >>> restrictions = [(1, 3), (4, 5)]

        IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
        (STARTS IN 0).
    """
    n_parts = min(len(mol1), len(mol2))
    mol1_ids_groups = _split_list(list(range(len(mol1))), n_parts)
    mol2_ids_groups = _split_list(list(range(len(mol2))), n_parts)
    restr = []
    for group1, group2 in zip(mol1_ids_groups, mol2_ids_groups):
        restr += [(i+offset1, j+offset2) for i in group1 for j in group2]
    return restr


def _split_list(alist, wanted_parts):
    """
    Taken from:
    "https://stackoverflow.com/questions/752308/split-list-into-smaller-lists"

    """
    length = len(alist)
    return [alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
            for i in range(wanted_parts)]
