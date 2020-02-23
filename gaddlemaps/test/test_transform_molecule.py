#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Test for the features from _transform_molecule file.

'''

import pytest
import numpy as np

from gaddlemaps import find_atom_random_displ, move_mol_atom


@pytest.fixture
def atom_pos():
    """Atoms positions
    """
    atom = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [2, 0, 0],
        [2, 1, 0],
        [3, 0, 0],
    ])
    return atom


@pytest.fixture
def bonds_info():
    """Bonds information
    """
    bonds = {
        0: [
            (1, 1.),
        ],
        1: [
            (0, 1.),
            (2, 1.),
            (3, 2**.5),
        ],
        2: [
            (1, 1.),
            (4, 1.),
        ],
        3: [
            (1, 2**.5),
        ],
        4: [
            (2, 1.),
        ],
    }

    return bonds


def test_find_atom_random_displ(atom_pos, bonds_info):
    """
    Test the direction of the found displacements.
    """
    # Normal to OX
    one_bond = find_atom_random_displ(atom_pos, bonds_info, 0)
    # Normal to OX
    two_bond = find_atom_random_displ(atom_pos, bonds_info, 2)
    # Normal to OX and OY
    three_bond = find_atom_random_displ(atom_pos, bonds_info, 1)
    assert not np.dot(one_bond, np.array([1, 0, 0]))
    assert not np.dot(two_bond, np.array([1, 0, 0]))
    assert not np.dot(three_bond, np.array([1, 0, 0]))
    assert not np.dot(three_bond, np.array([0, 1, 0]))


def test_move_mol_atom(atom_pos, bonds_info):
    """
    Test deformation of molecule in x dir
    """
    disp = np.array([1, 0, 0])
    new = move_mol_atom(atom_pos, bonds_info, 4, disp)
    new_guess = np.array([
        [1, 0, 0],
        [2, 0, 0],
        [3, 0, 0],
        [2, 1, 0],
        [4, 0, 0],
    ])
    assert np.array_equal(new, new_guess)
