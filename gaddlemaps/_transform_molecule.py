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
This module contains functions that change the molecules conformation.
'''


from collections import deque
from typing import Dict, List, Tuple, Deque

import numpy as np


def move_mol_atom(atoms_pos: np.ndarray,
                  bonds_info: Dict[int, List[Tuple[int, float]]],
                  atom_index: int = None,
                  displ: np.ndarray = None,
                  sigma_scale: float = 0.5) -> np.ndarray:
    """
    Moves an atom of a molecule respecting almost all bond distances.

    By default, a random atom is picked from atom_pos and moved randomly in a
    certain direction (see the published article for a better description of
    this step).

    Parameters
    ----------
    atoms_pos : numpy.ndarray
        An array with the positions of the molecule atoms in rows.
    bonds_info : dictionary
        A dict with the information of the bonds. Example:
            bonds_info = {0:[(1, 5.4), (2, 6.4)], 1:[(0, 5.4), ], ...}
        The keys refers to atom index and values are lists with tuples. Each
        tuple contains the bonded atom index and the bond length.
    atom_index : integer (Optional)
        The index of the atom to move (respecting the index of atoms_pos). If
        None is given a random one is taken.
    displ : numpy.ndarray (Optional)
        The displacement vector. If None is given a random displacement is
        calculated in the normal plane to the line jointing most nearest atoms.
    sigma_scale : float
        A factor to scale the sigma of the distribution of the
        displacement module.

    Returns
    -------
    modified_atoms_pos : numpy.ndarray
        An array with the modified positions of the atoms in rows.

    """

    atoms_pos = np.copy(atoms_pos)
    n_atoms = len(atoms_pos)
    if atom_index is None:
        atom_index = np.random.randint(n_atoms)

    if displ is None:
        displ = find_atom_random_displ(atoms_pos, bonds_info, atom_index,
                                       sigma_scale=sigma_scale)
    # Create a queue with the atoms to move to restore the bonds
    wait_queue = deque(range(n_atoms))
    # Remove an atom index if it was already moved
    atoms_pos[atom_index] += displ
    wait_queue.remove(atom_index)

    queue: Deque[Tuple[int, int, int]] = deque()
    for i in bonds_info[atom_index]:
        queue.append((atom_index, i[0], i[1]))  # type: ignore
        wait_queue.remove(i[0])

    while queue:
        ind1, ind2, bond = queue.pop()
        diferencia = atoms_pos[ind1] - atoms_pos[ind2]
        modulo = np.linalg.norm(diferencia)
        unit = diferencia/modulo
        atoms_pos[ind2] = atoms_pos[ind2] + (modulo - bond) * unit
        for bonds in bonds_info[ind2]:
            if bonds[0] in wait_queue:
                queue.append((ind2, bonds[0], bonds[1]))  # type: ignore
                wait_queue.remove(bonds[0])
    return atoms_pos


def find_atom_random_displ(atoms_pos: np.ndarray,
                           bonds_info: Dict[int, List[Tuple[int, float]]],
                           atom_index: int,
                           sigma_scale: float = 0.5) -> np.ndarray:
    """
    Finds a random displacement for the atom with a given index.

    This displacement is chosen in a perpendicular direction according to the
    number of bonded atoms.

    Parameters
    ----------
    atoms_pos : numpy.ndarray
        An array with the positions of the molecule atoms in rows.
    bonds_info : dictionary
        A dict with the information of the bonds. Example:
            bonds_info = {0:[(1, 5.4), (2, 6.4)], 1:[(0, 5.4), ], ...}
        The keys refers to atom index and values are lists with tuples. Each
        tuple contains the bonded atom index and the bond length.
    atom_index : integer (Optional)
        The index of the atom to move (respecting the index of atoms_pos). If
        None is given a random one is taken.
    sigma_scale : float
        A factor to scale the sigma of the distribution of the displacement
        module.

    Returns
    -------
    displ : numpy.ndarray
        The displacement vector to sum to the position of the interest atom.

    """
    n_bonded_ref = len(bonds_info[atom_index])
    # Calc the width of the displacements distribution
    sigma = bonds_info[atom_index][0][1] * sigma_scale

    if n_bonded_ref == 1:
        direction = np.cross(np.random.rand(3),
                             atoms_pos[bonds_info[atom_index][0][0]] -
                             atoms_pos[atom_index])
    elif n_bonded_ref == 2:
        direction = np.cross(np.random.rand(3),
                             atoms_pos[bonds_info[atom_index][0][0]] -
                             atoms_pos[bonds_info[atom_index][1][0]])
    elif n_bonded_ref >= 3:
        direction = np.cross(atoms_pos[bonds_info[atom_index][0][0]] -
                             atoms_pos[bonds_info[atom_index][2][0]],
                             atoms_pos[bonds_info[atom_index][0][0]] -
                             atoms_pos[bonds_info[atom_index][1][0]])
        direction *= np.random.choice([-1, 1])
    direction = direction / np.linalg.norm(direction)
    displ = direction * np.random.normal(0, sigma)
    return displ
