# -*- coding: utf-8 -*-
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
'''
This module provides the functor that perform the molecule transformation
in GADDLE MAPS.
'''

from collections import defaultdict
from typing import DefaultDict, Dict, List, Tuple

import numpy as np
from scipy.spatial.distance import euclidean

from ._auxilliary import calcule_base
from .components import Atom, Molecule


class ExchangeMap:
    """
    Functor to extrapolate atomic resolution to other configurations.

    When this class is initialized, a functor is created. It has to be
    initialized with the molecules in the initial and final resolution
    overlapped. Then you can call this method with new molecules in the
    initial resolution to obtain its representation in the final one. The new
    coordinates of the extrapolated molecules are scaled by the "scale_factor".
    This factor should be smaller than one if you want to extrapolate a complete
    system for future simulations. This avoids molecular overlapping and
    prevents the simulations to crash.

    Parameters
    ----------
    refmolecule : Molecule
        Molecule in the initial resolution.
    targetmolecule : Molecule
        Molecule in the final resolution.
    scale_factor : float, Optional
        The factor that modulates the scale of the atoms positions in the new
        resolution respect it closest atom in the initial resolution.

    Examples
    --------
       >>> BmimCG_ref = Molecule(fgro, fitp)
       >>> BmimAA_target = Molecule(groAA, itpAA)

       >>> transformation = ExchangeMap(BmimCG_ref, BmimAA_target)

       >>> BmimCG_new = Molecule(fgro2, fitp)
       >>> BmimAA_new = transformation(BmimCG_new)

    """

    def __init__(self, refmolecule: Molecule, targetmolecule: Molecule,
                 scale_factor: float = 0.5):
        self._refmolecule = refmolecule
        self._targetmolecule = targetmolecule
        self.scale_factor = scale_factor
        self._refsystems: Dict[int, Tuple[Tuple[np.ndarray, ...],
                                          np.ndarray]] = {}
        self._equivalences: Dict[int, int] = {}
        self._target_coordinates: Dict[int, np.ndarray] = {}
        self._calculate_refsystems(self._refmolecule)
        self._make_map()

    def _calculate_refsystems(self, molecule: Molecule):
        """
        Calculate the referencce system of all the atoms that will
        be used as reference points for the coordiante
        transformation

        """
        n_atoms = len(molecule)
        # Special cases
        if n_atoms in [1, 2]:
            pos = molecule.atoms_positions
            rand_pos = [np.random.rand(3) + pos[0] for _ in range(3-n_atoms)]
            positions = np.append(pos, rand_pos, axis=0)
            coord_syst = calcule_base(positions)
            self._refsystems[hash(molecule[0])] = coord_syst
        else:
            self._calculate_refsystems_general(molecule)

    def _calculate_refsystems_general(self, molecule: Molecule):
        """
        Version of _calculate_refsystems for molecules of 3 or more atoms.

        """
        # Select the atoms that are bonded at least to 2 other atoms
        for atom in molecule:
            if len(atom.bonds) >= 2:
                # Get 2 of the atoms that are bonded to the atom
                ind1, ind2 = atom.closest_atoms()
                neighbour1, neighbour2 = (molecule[ind1], molecule[ind2])
                positions = [atom.position, neighbour1.position,
                             neighbour2.position]
                # se the 3 atoms to calculate the base
                coord_syst = calcule_base(positions)
                self._refsystems[hash(atom)] = coord_syst

    def _make_map(self):
        """
        Make the conections between the atoms of the target and the references.

        It also calculates the coordinates of target in the ref coordinates
        system.

        """
        for atom in self._targetmolecule:
            name = self._find_closest_ref(atom)
            proyection = self._proyect_point(name, atom)
            self._equivalences[hash(atom)] = name
            self._target_coordinates[hash(atom)] = proyection

    def _find_closest_ref(self, targetatom: Atom) -> int:
        """
        Returns the name of the closest reference system to the atom

        """
        def ref_pos(index):
            return self._refmolecule[index].position

        distances = [(euclidean(targetatom.position, ref_pos(index)), index)
                     for index in self._refsystems]
        return sorted(distances)[0][1]

    def _restore_point(self, atomref: int,
                       proyection: np.ndarray) -> np.ndarray:
        """
        Given an atomref and a projection restores the coordinates of the
        projected point.

        """
        center = self._refsystems[atomref][1]
        vectores = np.array(self._refsystems[atomref][0])
        return center + np.dot(proyection, vectores)

    def _proyect_point(self, atomref: int, atomtarget: Atom) -> np.ndarray:
        origen = self._refsystems[atomref][1]
        proyect = atomtarget.position - origen
        proyect = np.dot(self._refsystems[atomref][0], proyect)
        return proyect * self.scale_factor

    def _restore_molecule(self) -> Molecule:
        new_mol = self._targetmolecule.copy()
        for atom in new_mol:
            refname = self._equivalences[hash(atom)]
            proyection = self._target_coordinates[hash(atom)]
            atom.position = self._restore_point(refname, proyection)
        return new_mol

    def __call__(self, refmolecule: Molecule) -> Molecule:
        """
        This function takes as argument a molecule like the refmolecule, but
        in other position and returns its targetmolecule equivalent as a new
        molecule with the same res_number than the input.

        """
        if not isinstance(refmolecule, Molecule):
            raise TypeError("Argument must be a Molecule")
        if self._refmolecule != refmolecule:
            raise TypeError(("refmolecule must be:\n{}"
                             "").format(self._refmolecule))

        self._calculate_refsystems(refmolecule)
        new_mol = self._restore_molecule()
        # change the indexes
        new_mol.resids = refmolecule.resids
        return new_mol

    @property
    def equivalences(self) -> Dict[int, List[int]]:
        """
        dict of int to list of int : {r1_atom_index: [closest_r2_atoms_indexs]}
        """
        rev:  DefaultDict[int, List[int]] = defaultdict(list)
        for ind_r1, ind_r2 in self._equivalences.items():
            rev[ind_r2].append(ind_r1)
        return dict(rev)
