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
Test for the functions defined in the _backend submodule.
'''

import os

import pytest
import numpy as np

from gaddlemaps.components import System, Molecule
from gaddlemaps import (Chi2Calculator, accept_metropolis,
                        check_backend_installed)
from gaddlemaps._backend import _minimize_molecules


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def molecule_aa() -> Molecule:
    """
    Molecule instance of curcumine all atom.
    """
    fgro = os.path.join(ACTUAL_PATH, '../data/CUR_AA.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/CUR_AA.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def molecule_cg() -> Molecule:
    """
    Molecule instance of curcumine coarse-grained.
    """
    fgro = os.path.join(ACTUAL_PATH, '../data/CUR_map.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/CUR_CG.itp')
    return System(fgro, fitp)[0]


def test_chi2_molecules():
    """
    Test chi2 calculations with simple molecules.
    """
    mol1 = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
    ])

    mol2 = np.array([
        [0, 0, 1],
        [0, 2, 0],
    ])

    chi2_test = 1 + 2 + 1
    chi2_molecules = Chi2Calculator(mol1, mol2)
    assert chi2_test == chi2_molecules(mol2)
    # With restrictions one atom does not take role
    restr = np.array([(0, 0), (1, 0), (2, 0)])
    mol2 = np.array([
        [0, 0, 1],
        [0, 9, 0],
    ])
    chi2_test = 1 + 2 + 2
    chi2_test *= 1.1
    chi2_molecules = Chi2Calculator(mol1, mol2, restr)
    assert chi2_test == chi2_molecules(mol2)


def test_chi2_molecules_inv():
    """
    Same test as test_chi2_molecules but changing the order in input
    molecules.
    """
    mol1 = np.array([
        [0, 0, 1],
        [0, 2, 0],
    ])

    mol2 = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
    ])
    restr = np.array([(0, 0), (1, 1), (2, 1)])
    chi2_molecules = Chi2Calculator(mol2, mol1, restr)
    assert chi2_molecules(mol1) == 7


def test_minimize_molecules(molecule_aa: Molecule,
                            molecule_cg: Molecule):
    """
    Test the molecule minimization backend
    """
    molecules = [molecule_aa, molecule_cg]
    mol1_positions = molecules[0].atoms_positions
    mol2_positions = molecules[1].atoms_positions
    mol2_bonds_info = molecules[1].bonds_distance
    mol2_com = molecules[1].geometric_center
    n_steps = 5000*len(molecules[1])
    mol1_bonds_info = molecules[0].bonds_distance.values()
    translation_width = 2*min(t[1] for l in mol1_bonds_info for t in l)
    deformation_types = (0, 1, 2)
    mol2_positions_new = _minimize_molecules(mol1_positions, mol2_positions,
                                             mol2_com, .5,
                                             n_steps, [(0, 0)],
                                             mol2_bonds_info, translation_width,
                                             deformation_types)
    assert mol2_positions_new.shape == mol2_positions.shape
    if check_backend_installed():
        from cython_backend._backend import py_minimize_molecules
        positions = py_minimize_molecules(mol1_positions, mol2_positions,
                                          mol2_com, .5,
                                          n_steps, [(0, 0)],
                                          mol2_bonds_info, translation_width,
                                          deformation_types)
        mol2_positions_new_cpp = np.array(positions)
        assert np.all(np.isclose(mol2_positions_new, mol2_positions_new_cpp,
                                 rtol=1, atol=1))


def test_accept_metropolis():
    """
    Test the acceptance algorithm with non random cases.
    """
    assert accept_metropolis(50, 40)
    assert accept_metropolis(40, 40)
    assert not accept_metropolis(30, 40, acceptance=0)
    assert not accept_metropolis(40, 50)
