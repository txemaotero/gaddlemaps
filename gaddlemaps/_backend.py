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
"""
This module implements a pure python version of the minimization
algorithm to use if the c++ one is not available.
"""

import sys
from typing import Dict, List, Optional, Tuple
import warnings

import numpy as np
from scipy.spatial.distance import cdist

from . import move_mol_atom
from ._auxilliary import rotation_matrix


def check_backend_installed(warn_missing=False) -> bool:
    """
    Returns wether the compiled backend is installed and accesible or not. It
    can also raisea warning if the backend is not fond.

    Parameters
    ----------
    warn_missing: Optional[bool]
        If True the function will also raise a warning if the backend is not
        installed. If False no warnings will be raised. Default False

    Returns
    -------
    installed: bool
        Whether the compiled backend is installed and accesible or not.

    """
    try:
        from cython_backend._backend import py_minimize_molecules
        return True
    except ImportError:
        if warn_missing:
            text = ("Compiled backend was not found, GADDLE MAPS will fallback"
                    " to the python aligmnet engine. You will experiencie "
                    "degraded performance in the aligment, but the results will"
                    " be the same. If the current performance is good enough"
                    " you may ignore this warning. In any case, refer to the "
                    "documentation for the installation guide for "
                    "the compiled backed.")
            warnings.warn(text)
        return False


def minimize_molecules(mol1_positions: np.ndarray,
                       mol2_positions: np.ndarray,
                       mol2_com: np.ndarray,
                       sigma_scale: float,
                       n_steps: int,
                       restriction: List[Tuple[int, int]],
                       mol2_bonds_info: Dict[int, List[Tuple[int, float]]],
                       displacement_module: float,
                       sim_type: Tuple[int, ...]):
    """
    Minimizes the distance between two molecules.

    Parameters
    ----------
    mol1_positions : np.array((N, 3))
        The positions of the atoms of the static molecule.
    mol2_positions : np.array((N, 3))
        The positions of the atoms of the mobile molecule.
    mol2_com : np.array(3)
        The center of mass (or the geometric center) of the mobile molecule.
        This will be used to apply de rotations.
    sigma_scale : float
        A number that module the amplitude of the single atom displacements.
    n_steps : int
        The number of steps without changes in the chi2 in the Monte-Carlo
        alignment.
    restriction : list of tuple of int
        A list of tuples with pairs of atom numbers corresponding to mol1
        and mol2 atoms molecules. The align will be performed privileging
        configurations where those atoms are close.

        Example
        -------
        >>> restrictions = [(1, 3), (4, 5)]

        **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
        (STARTS IN 0).**

    mol2_bonds_info : defaultdict of int: list of tuple(int, float)
        A complex data structure that collect the information of the bond
        distances. The key of the property corresponds to the atom index in the
        molecule. The value is a list with tuples. For each tuple, the first
        value corresponds with the index of the bonded atom and the second is
        the length of the bond.
    same_com : bool
        If True, translations are not allowed.
    anchura : float
        A number that modules the amplitude in the translations displacements.
    sim_type : tuple of int
        Specifies the type of the minimization. Possible options:
             0 : Translation
             1 : Rotation
             2 : Individual atom move
        If it is a list, all specified methods are combined.

    """
    if check_backend_installed(warn_missing=True):
        from cython_backend._backend import py_minimize_molecules
        positions = py_minimize_molecules(mol1_positions, mol2_positions,
                                          mol2_com, sigma_scale, n_steps,
                                          restriction, mol2_bonds_info, displacement_module,
                                          sim_type)
        return np.array(positions)

    mol2_positions = _minimize_molecules(mol1_positions, mol2_positions,
                                         mol2_com, sigma_scale, n_steps,
                                         restriction, mol2_bonds_info,
                                         displacement_module, sim_type)
    return mol2_positions


def _minimize_molecules(mol1_positions: np.ndarray,
                        mol2_positions: np.ndarray,
                        mol2_com: np.ndarray,
                        sigma_scale: float,
                        n_steps: int,
                        restriction: List[Tuple[int, int]],
                        mol2_bonds_info: Dict[int, List[Tuple[int, float]]],
                        displacement_module: float,
                        sim_type: Tuple[int, ...]):
    # Performance stuffs
    _chi2_molecules = Chi2Calculator(mol1_positions, mol2_positions,
                                     restriction)
    _move_mol_atom = move_mol_atom
    _accept_metropolis = accept_metropolis
    _mean = np.mean
    _rand_unif = np.random.uniform
    _rand_norm = np.random.normal
    _rotation_matrix = rotation_matrix
    _choice = np.random.choice
    _dot = np.dot

    # Init chi2
    chi2 = _chi2_molecules(mol2_positions)
    chi2_min = chi2
    counter = 0

    while counter < n_steps:
        change = _choice(sim_type)

        # Translations
        if change == 0:
            desplazamiento = _rand_norm(0, displacement_module, 3)
            test = mol2_positions + desplazamiento
            chi2_new = _chi2_molecules(test)

        # Rotations: random with angles between 0 y pi
        elif change == 1:
            # geometric center
            mol2_com = _mean(mol2_positions, axis=0)
            axis = _rand_unif(-1, 1, 3)
            theta = _rand_norm(0, np.pi/4.)
            rot_matrix = _rotation_matrix(axis, theta)
            # test configuration
            test = _dot(mol2_positions-mol2_com, rot_matrix) + mol2_com
            chi2_new = _chi2_molecules(test)

        # Single atom displacement
        elif change == 2:
            test = _move_mol_atom(mol2_positions, mol2_bonds_info,
                                  sigma_scale=sigma_scale)
            chi2_new = _chi2_molecules(test)

        # Check if the new config has to be accepted
        if _accept_metropolis(chi2, chi2_new):
            mol2_positions = test
            chi2 = chi2_new

            # Check if the counter has to be restarted
            if chi2 < chi2_min:
                chi2_min = chi2
                sys.stdout.write('\r')
                sys.stdout.write('\tChi2 = %10.9f' % (chi2_min))
                sys.stdout.flush()
                counter = 0
                continue
        counter += 1

    print('\n')
    return mol2_positions


def accept_metropolis(energy_0: float, energy_1: float,
                      acceptance: float = 0.01):
    """
    Evaluates if the new configuration with energy_1 is accepted or not.

    This function returns True if the configuration with energy_1 is accepted
    and False otherwise. If energy_1 <= energy_0 the configuration is accepted.
    Else, it will be also accepted if a random number between 0 and 1 is
    smaller than acceptance*energy_0/energy_1.

    Parameters
    ----------
    energy_0 : float
        The value of the energy before the change.
    energy_1 : float
        The value of the energy after the change.
    acceptance : float (Optional)
        acceptance factor that regulates how many unfavorable cases are accepted

    Returns
    -------
    change : Bool
        True if the change is accepted and False if not.

    """

    factor = energy_0/energy_1
    condition = factor >= 1
    if condition:
        return condition
    return np.random.rand() <= acceptance*factor


class Chi2Calculator:
    """
    Functor that calculates the chi2 between 2 positions arrays.

    This class is initialized with the restrictions that must be taken into
    account. This avoid repeating array accessing. Once the object is
    initialized it can be called with new atomic positions  to calculate the new
    chi2 value associated with the distance between set of points. This distance
    is calculated quite differently depending on the given restrictions. In case
    of no restrictions see "chi2_molecules" method. If some restrictions are
    given the distances between atoms is calculated between of pairs of atoms
    specified in the restrictions instead of between closest atoms.

    Parameters
    ----------
    mol1 : numpy.ndarray((N, 3))
        Array with the atomic positions of the molecule that will be still.
        These positions will remain constant in  future object call.
    mol2 : numpy.ndarray((N, 3))
        Array with the atomic positions of the molecule that will be changing
        during the optimum overlap finding process.
    restriction : numpy.ndarray((N, 2)) or array convertible, optional
        A list of tuples with pairs of atom numbers corresponding to mol1
        and mol2 atoms molecules. The align will be performed privileging
        configurations where those atoms are close. By default is set to [].

        Example
        -------
        >>> restrictions = [(1, 3), (4, 5)]

        **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE**

    """

    def __init__(self, mol1: np.ndarray, mol2: np.ndarray,
                 restrictions: Optional[np.ndarray] = None):
        self._mol1_positions = mol1
        if restrictions is None or len(restrictions) == 0:
            self._meth_to_call = self.chi2_molecules
        else:
            self.restrictions = np.array(restrictions)

            restriction1, self.restriction2 = self.restrictions.T
            self.set_restriction2 = set(self.restriction2)

            self._meth_to_call = self._chi2_molecules_with_restrains

            self.len_mol2 = len(mol2)

            mol1_not_restriction_mask = np.ones(len(mol1), dtype=np.bool)
            mol1_not_restriction_mask[restriction1] = False

            self._mol1_not_restriction = mol1[mol1_not_restriction_mask]
            self._mol1_restriction = mol1[restriction1]

            if not mol1_not_restriction_mask.any():
                self._meth_to_call = self._chi2_molecules_only_restrains
                self.n_cg_far_fact = 1.1**(self.len_mol2
                                           - len(self.set_restriction2))

    def __call__(self, mol2: np.ndarray) -> float:
        return self._meth_to_call(mol2)

    def _chi2_molecules_restrains_contrib(self, mol2: np.ndarray) -> float:
        """
        Calculates the contribution to Chi2 of the restrains.
        """
        mol2_restrictions = mol2[self.restriction2]
        return np.sum((self._mol1_restriction - mol2_restrictions)**2)

    def _chi2_molecules_only_restrains(self, mol2: np.ndarray) -> float:
        chi2 = self._chi2_molecules_restrains_contrib(mol2)
        return chi2 * self.n_cg_far_fact

    def _chi2_molecules_with_restrains(self, mol2: np.ndarray) -> float:
        """
        Computes the chi2 by calculating the distance between molecules.

        This function computes a chi2 based on the distance between nearest
        atoms of two molecules. Basically sums the distance between the atoms
        of mol2 and the nearest atom of mol1.

        Parameters
        ----------
        mol2 : numpy.ndarray((N,3))
            Position of the second molecule atoms.

        Returns
        -------
        chi2 : float
            Computed chi2.

        """

        # Easy distance calculation
        chi2 = self._chi2_molecules_restrains_contrib(mol2)
        # Distance  between all the atoms that are not in restrictions
        distances = cdist(self._mol1_not_restriction, mol2, 'sqeuclidean')
        chi2 += np.sum(distances.min(axis=1))

        # Penalize the value if some atom is not close to any other atom
        n_cg_far = (self.len_mol2 -
                    len(self.set_restriction2.union(distances.argmin(axis=1))))
        if n_cg_far:
            chi2 *= 1.1**n_cg_far
        return chi2

    def chi2_molecules(self, mol2: np.ndarray) -> float:
        """
        Computes the chi2 by calculating the distance between molecules.

        This function computes a chi2 based on the distance between nearest
        atoms of two molecules. Basically sums the distance between the atoms
        of mol2 and the nearest atom of mol1.

        Parameters
        ----------
        mol2 : numpy.ndarray((N,3))
            Position of the second molecule atoms.

        Returns
        -------
        chi2 : float
            Computed chi2.

        """

        distances = cdist(self._mol1_positions, mol2, 'sqeuclidean')
        chi2 = np.sum(distances.min(axis=1))
        n_cg_far = len(mol2) - len(set(distances.argmin(axis=1)))
        if n_cg_far:
            chi2 *= 1.1**n_cg_far
        return chi2
