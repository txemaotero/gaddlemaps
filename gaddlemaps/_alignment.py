# -*- coding: utf-8 -*-
'''
This module contains features to manage the alignment engine.
'''

import numpy as np

from typing import Tuple, Union, List, Optional, Any, Dict

from .components import Molecule
from .parsers import GroFile
from . import ExchangeMap, minimize_molecules


class Alignment(object):
    """
    Class to manage the molecules alignment.

    This class serve as interface during the molecular alignment process. It
    will take two molecules in different representations (e.g. coarse grained
    and atomistic), one of them will be taken as initial representation and the
    other as the final one. These molecules can be passed as arguments in the
    object initialization or assigned after its creation. 

    NOTE: Once both molecules are assigned you can only set them again with the
    same molecule type (although it may have different position).

    Once the molecules are fixed you can call the "aling_molecules" method to
    find the optimum overlap between both representation. Then, you should
    call the "init_exchange_map" method to generate the map to change between
    resolution for other start molecules configurations. You can also call
    the "write_comparative_gro" method to write a .gro file with the two
    aligned molecules with different residue names to visualize the found
    overlap in other visualization programs.

    Parameters
    ----------
    start : Molecule
        The molecule to align in the initial resolution.
    end : Molecule
        The molecule to align in the final resolution.

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

    def __init__(self, start: Optional[Molecule] = None,
                 end: Optional[Molecule] = None):
        super(Alignment, self).__init__()
        self._start = None
        self._end = None
        self.start = start
        self.end = end
        self.exchange_map: Optional[ExchangeMap] = None

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
    def start(self, molecule: Optional[Molecule]):
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
        Molecule : The molecule in the final resolution.

        Setter Raises
        -------------
        TypeError
            If the molecule is not instance of Molecule
        ValueError
            If the molecule does not match with the start.
        """
        return self._end

    @end.setter
    def end(self, molecule: Optional[Molecule]):
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

    def align_molecules(self, restrictions: Optional[List[Tuple[int, int]]] = None,
                        deformation_types: Optional[Tuple[int, ...]] = None,
                        ignore_hydrogens: bool = True,
                        auto_guess_protein_restrictions: bool = True):
        """
        Starts the alignment engine to find the optimal overlap between molecs.

        NOTE: If None is input as restrictions and the molecules to align
        have multiple residues and also the auto_guess_protein_restrictions
        parameter is True, they will be guessed by residue matching (see
        guess_protein_restrains function). However take into account that the
        time taken by this step scales with the square of the number of atoms if
        no restrictions are given and linearly if they are guessed for all the
        atoms by residue matching.

        If you know that the molecules in both resolution are in the same
        configuration you may want to perform the alignment just rotating and
        translating the molecules as a whole and avoid molecular deformation by
        setting the "deformation_types" parameter to (0, 1). This will translate
        in a better performance.

        On the other hand, in the most of the cases, the hydrogen atoms are not
        very relevant atoms for the alignment so they are by default (see
        "ignore_hydrogens" parameter) not taken into account in the distance
        minimization process when they belong to the molecule with larger number
        of atoms (the one that will remain still in the process).

        Parameters
        ----------
        restrictions : list of tuple of int, optional
            A list of tuples with pairs of atom numbers corresponding to
            start and end atoms indexes in the molecules. The align will be
            performed privileging configurations where those atoms are close.
            By default is set to [].

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
            is not moved in the alignment engine.
        auto_guess_protein_restrictions : bool, optional
            If True automatic restrictions will try to be guessed if the
            molecules to be aligned have multiple residues.

        Raises
        ------
        IOError
            If the automatic restrictions guess failed. In that case, consider
            to set "auto_guess_protein_restrictions" parameter to False.

        """

        if restrictions is None:
            # Check for multiple residues
            restrictions = []
            if len(self.start.resnames) > 1 and auto_guess_protein_restrictions:
                try:
                    restrictions = guess_protein_restrains(self.start,
                                                           self.end)
                except IOError:
                    raise IOError(('Automatic restrictions can not be guessed.'
                                   ' Try to set '
                                   'auto_guess_protein_restrictions to '
                                   'False.'))

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

    def init_exchange_map(self, scale_factor: float = 0.5):
        """
        Initializes the exchange map with the current molecules configuration.

        Parameters
        ----------
        scale_factor : float
            The compression factor to apply to mapped molecules (see the
            ExchangeMap class documentation for more information).

        """
        self.exchange_map = ExchangeMap(self.start, self.end, scale_factor)

    def write_comparative_gro(self, fname: Optional[str] = None):
        """
        Writes a .gro file with start and end molecules to check the overlap.

        The molecules have START and END residue names respectively to ease the
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
    def info(self) -> Dict[str, Any]:
        """
        dict : A dictionary with the needed information to restore the
            Alignment object.
        """
        ex_map = (None if self.exchange_map is None
                  else self.exchange_map.scale_factor)

        def mol_dic(mol): return mol if mol is None else mol.info
        info = {
            'start': mol_dic(self.start),
            'end': mol_dic(self.end),
            'exchange_map': ex_map,
        }
        return info

    @classmethod
    def from_info(cls, info_dict: Dict[str, Any]) -> 'Alignment':
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


def remove_hydrogens(molecule: Molecule,
                     restrictions: List[Tuple[int, int]]) -> Tuple[np.ndarray,
                                                                   List[Tuple[int, int]]]:
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


def guess_protein_restrains(mol1: Molecule,
                            mol2: Molecule) -> List[Tuple[int, int]]:
    """
    Guess restrains for molecules with multiple residues.

    Checks if mol1 and mol2 have the same residue names and create restrains
    pairing atoms with the same residues. This function only works with
    molecules with multiple residues.

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
        If the input molecules have different residue names.
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
    restrains: List[Tuple[int, int]] = []
    offset1 = 0
    offset2 = 0
    for mol_res1, mol_res2 in zip(mol1.molecules, mol2.molecules):
        restrains += guess_molecule_restrains(mol_res1, mol_res2, offset1,
                                              offset2)
        offset1 += len(mol_res1)
        offset2 += len(mol_res2)
    return restrains


def guess_molecule_restrains(mol1: Molecule, mol2: Molecule,
                             offset1: int = 0,
                             offset2: int = 0) -> List[Tuple[int, int]]:
    """
    Guess restrains for molecules with one residue.

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
    restr: List[Tuple[int, int]] = []
    for group1, group2 in zip(mol1_ids_groups, mol2_ids_groups):
        restr += [(i+offset1, j+offset2) for i in group1 for j in group2]
    return restr


def _split_list(alist: List[Any], wanted_parts: int) -> List[List[Any]]:
    length = len(alist)
    return [alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
            for i in range(wanted_parts)]
