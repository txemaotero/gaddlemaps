'''
Test for the functions defined in the _backend submodule.
'''

import numpy as np

from gaddlemaps import Chi2Calculator, accept_metropolis


def test_chi2_molecules():
    """
    Test simple cases with simple molecules.
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


def test_accept_metropolis():
    """
    Test the deterministic case.
    """
    assert accept_metropolis(50, 40)
    assert accept_metropolis(40, 40)
    assert not accept_metropolis(30, 40, acceptance=0)
    assert not accept_metropolis(40, 50)
