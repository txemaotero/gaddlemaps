#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
This moldule contains functions that change the molecules conformation.
'''


from collections import deque

import numpy as np


def move_mol_atom(atoms_pos, bonds_info, atom_index=None, displ=None,
                  sigma_scale=0.5):
    """
    This function moves an atom of a molecule respecting almost all bond
    distances.

    Parameters
    ----------
    atoms_pos : numpy.ndarray
        An array with the positions of the molecule atoms in rows.
    bonds_info : dictionary
        A dict with the information of the bonds. Example:
            bonds_info = {0:[(1, 5.4), (2, 6.4)], 1:[(0, 5.4), ], ...}
        The keys refers to atom index and values are lists with tupples. Each
        tupple contains the bonded atom index and the bond length.
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
        An array with the modified positions of the molecule atoms in rows.

    """

    # Se crea una copia de las posiciones
    atoms_pos = np.copy(atoms_pos)
    # Se mira si hay que coger in random
    n_atoms = len(atoms_pos)
    if atom_index is None:
        atom_index = np.random.randint(n_atoms)
    # Se encuentra el desplazamiento aleatorio de ser necesario
    if displ is None:
        displ = find_atom_random_displ(atoms_pos, bonds_info, atom_index,
                                       sigma_scale=sigma_scale)
    # Comienza el movimiento de la molécula Se crea una lista de espera con los
    # índices de los átomos que aún no se han movido
    espera = deque(range(n_atoms))
    # Se mueve el primero y se elimina de la lista
    atoms_pos[atom_index] += displ
    espera.remove(atom_index)
    # Se inicia la cola que contendrá los índices de los átomos que se van a
    # mover y en qué orden
    cola = deque()
    # Se recorren los átomos enlazados al primero y se añaden a la cola y se
    # eliminan de espera
    for i in bonds_info[atom_index]:
        cola.append((atom_index, i[0], i[1]))
        espera.remove(i[0])
    # Se van moviendo todos
    while cola:
        # Se calcula la nueva posición de átomo que le toque
        ind1, ind2, bond = cola.pop()
        diferencia = atoms_pos[ind1] - atoms_pos[ind2]
        modulo = np.linalg.norm(diferencia)
        unit = diferencia/modulo
        atoms_pos[ind2] = atoms_pos[ind2] + (modulo - bond) * unit
        # Se añaden a la cola los átomos enlazados que no se han considerado
        for bonds in bonds_info[ind2]:
            if bonds[0] in espera:
                cola.append((ind2, bonds[0], bonds[1]))
                espera.remove(bonds[0])
    return atoms_pos


def find_atom_random_displ(atoms_pos, bonds_info, atom_index, sigma_scale=0.5):
    """
    This function finds a random displacement of the atom with the given index
    in a perpendicular direction according to its bonded atom.

    Parameters
    ----------
    atoms_pos : numpy.ndarray
        An array with the positions of the molecule atoms in rows.
    bonds_info : dictionary
        A dict with the information of the bonds. Example:
            bonds_info = {0:[(1, 5.4), (2, 6.4)], 1:[(0, 5.4), ], ...}
        The keys refers to atom index and values are lists with tupples. Each
        tupple contains the bonded atom index and the bond length.
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
    # Se lee la anchura de la distr de desplazamientos
    sigma = bonds_info[atom_index][0][1] * sigma_scale
    # Se analizan los casos posibles
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
