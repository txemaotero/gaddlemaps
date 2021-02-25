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
Test for the features from _transform_molecule submodule.
'''

from typing import Dict, List, Tuple

import pytest
import numpy as np

from gaddlemaps import find_atom_random_displ, move_mol_atom


@pytest.fixture
def atom_pos() -> np.ndarray:
    """
    Some artificial atom positions
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
def bonds_info() -> Dict[int, List[Tuple[int, float]]]:
    """
    Some artificial bonds information
    """
    bonds = {
        0: [ (1, 1.), ],
        1: [ (0, 1.), (2, 1.), (3, 2**.5), ],
        2: [ (1, 1.), (4, 1.), ],
        3: [ (1, 2**.5), ],
        4: [ (2, 1.), ],
    }
    return bonds


def test_find_atom_random_displ(atom_pos: np.ndarray,
                                bonds_info: Dict[int, List[Tuple[int, float]]]):
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


def test_move_mol_atom(atom_pos: np.ndarray,
                       bonds_info: Dict[int, List[Tuple[int, float]]]):
    """
    Test deformation of molecule in the x direction.
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
