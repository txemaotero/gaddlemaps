# -*- coding: utf-8 -*-
#!/usr/bin/python

"""
This module implements a pure python version of the minimization
algorithm to use if the c++ one is not available.
"""


from __future__ import division


import sys
import numpy as np
from scipy.spatial.distance import cdist
from . import move_mol_atom
from ._auxilliary import rotation_matrix


def minimize_molecules(mol1_positions, mol2_positions, mol2_com, sigma_scale,
                       n_steps, restriction, mol2_bonds_info, anchura,
                       sim_type):
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

        Example:
        >>> restrictions = [(1, 3), (4, 5)]

        IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
        (STARTS IN 0).

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
    not_repeat_calculus : bool
        Unused input.

    """
    try:
        from cython_backend._backend import py_minimize_molecules
        positions = py_minimize_molecules(mol1_positions, mol2_positions,
                                          mol2_com, sigma_scale, n_steps,
                                          restriction, mol2_bonds_info, anchura,
                                          sim_type)
        return np.array(positions)
    except ImportError:
        pass

    # Performance stufs
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
        # Se le escoge que se le hace a la molecula
        cambio = _choice(sim_type)

        # Translations
        if cambio == 0:
            desplazamiento = _rand_norm(0, anchura, 3)
            test = mol2_positions + desplazamiento
            chi2_new = _chi2_molecules(test)

        # Rotations: Estas son aleatorias con ángulos entre 0 y pi
        elif cambio == 1:
            # Se calcula el centro geometrico
            mol2_com = _mean(mol2_positions, axis=0)
            axis = _rand_unif(-1, 1, 3)
            theta = _rand_norm(0, np.pi/4.)
            # Se obtiene la matriz de rotación
            matriz_rot = _rotation_matrix(axis, theta)
            # Se rota
            test = _dot(mol2_positions-mol2_com, matriz_rot) + mol2_com
            chi2_new = _chi2_molecules(test)

        # Single atom displacement
        elif cambio == 2:
            test = _move_mol_atom(mol2_positions, mol2_bonds_info,
                                  sigma_scale=sigma_scale)
            chi2_new = _chi2_molecules(test)

        # Se mira si se acepta la conf.
        if _accept_metropolis(chi2, chi2_new):
            mol2_positions = test
            chi2 = chi2_new

            # Se mira si hay que reiniciar el counter
            if chi2 < chi2_min:
                chi2_min = chi2
                sys.stdout.write('\r')
                sys.stdout.write('\tChi2 = %10.9f'%(chi2_min))
                sys.stdout.flush()
                counter = 0
                continue
        counter += 1

    print('\n')
    return mol2_positions


def accept_metropolis(energy_0, energy_1, acceptance=0.01):
    """
    Evaluates the metropolis algorithm to check if one configuration is
    accepted or not.

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


class Chi2Calculator(object):
    """
    Functor that calculates the chi2 between 2 positions arrays.

    This class is intitialized with the restrictions that must be taken into
    account. This avoid repeating array accessing.

    Parameters
    ----------
    restriction : numpy.ndarray((N, 2)) or array convertible, optional
        A list of tuples with pairs of atom numbers corresponding to mol1
        and mol2 atoms molecules. The align will be performed privileging
        configurations where those atoms are close. By default is set to [].

        Example:
        >>> restrictions = [(1, 3), (4, 5)]

        IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE

    """
    def __init__(self, mol1, mol2, restrictions=None):
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

    def __call__(self, mol2):
        return self._meth_to_call(mol2)

    def _chi2_molecules_restrains_contrib(self, mol2):
        """
        Calculates the contribution to Chi2 of the restrains.
        """
        mol2_restrictions = mol2[self.restriction2]
        return np.sum((self._mol1_restriction - mol2_restrictions)**2)

    def _chi2_molecules_only_restrains(self, mol2):
        chi2 = self._chi2_molecules_restrains_contrib(mol2)
        return chi2 * self.n_cg_far_fact

    def _chi2_molecules_with_restrains(self, mol2):
        """
        Computes the chi2 by calculating the distance between molecules.

        This function computes a chi2 based on the distance between nearest
        atoms of two molecules. Basically sums the distance between the atoms
        of mol2 and the nearest atom of mol1.

        Parameters
        ----------
        mol1 : numpy.ndarray((N,3))
            Position of the first molecule atoms.
        mol2 : numpy.ndarray((N,3))
            Position of the second molecule atoms.
        restriction : numpy.ndarray((N, 2)), optional
            A list of tuples with pairs of atom numbers corresponding to mol1
            and mol2 atoms molecules. The align will be performed privileging
            configurations where those atoms are close. By default is set to [].

            Example:
            >>> restrictions = [(1, 3), (4, 5)]

            IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE

        Returns
        -------
        chi2 : float
            Computed chi2.

        """

        chi2 = self._chi2_molecules_restrains_contrib(mol2)

        # Se calculan las distancias entre pares de átomos
        distances = cdist(self._mol1_not_restriction, mol2, 'sqeuclidean')

        # hay que sumar la distancia mínima de cada fila
        chi2 += np.sum(distances.min(axis=1))

        # Se sube su valor si hay un CG que no está cerca de ninguno
        n_cg_far = (self.len_mol2 -
                    len(self.set_restriction2.union(distances.argmin(axis=1))))
        if n_cg_far:
            chi2 *= 1.1**n_cg_far
        return chi2

    def chi2_molecules(self, mol2):
        """
        Computes the chi2 by calculating the distance between molecules.

        This function computes a chi2 based on the distance between nearest
        atoms of two molecules. Basically sums the distance between the atoms
        of mol2 and the nearest atom of mol1.

        Parameters
        ----------
        mol1 : numpy.ndarray((N,3))
            Position of the first molecule atoms.
        mol2 : numpy.ndarray((N,3))
            Position of the second molecule atoms.

        Returns
        -------
        chi2 : float
            Computed chi2.

        """

        # Se calculan las distancias entre pares de átomos
        distances = cdist(self._mol1_positions, mol2, 'sqeuclidean')
        # hay que sumar la distancia mínima de cada fila
        chi2 = np.sum(distances.min(axis=1))
        # Se sube su valor si hay un CG que no está cerca de ninguno
        n_cg_far = len(mol2) - len(set(distances.argmin(axis=1)))
        if n_cg_far:
            chi2 *= 1.1**n_cg_far
        return chi2
