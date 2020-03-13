"""
Tests for the _top_parsers submodule.
"""

import os

import pytest

from gaddlemaps.parsers import ItpFile, itp_top
from gaddlemaps.parsers._top_parsers import (_itp_top_atoms, _itp_top_name,
                                             _parse_itp_bonds)

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def bf4_itp_fname() -> str:
    """
    File name with the itp of a BF4 molecule.
    """
    return os.path.join(ACTUAL_PATH, '../../data/BF4_AA.itp')


@pytest.fixture
def bf4_itp_file(bf4_itp_fname: str) -> ItpFile:
    """
    Loaded ItpFile object with a BF4 molecule.
    """
    return ItpFile(bf4_itp_fname)


def test_itp_top_name(bf4_itp_file: ItpFile):
    """
    Test for the _itp_top_name function.
    """
    assert _itp_top_name(bf4_itp_file) == 'BF4'
    

def test_parse_itp_bonds(bf4_itp_file: ItpFile):
    """
    Test for the _parse_itp_bonds function.
    """
    bonds = _parse_itp_bonds(bf4_itp_file)
    test_bonds = [(1, 2), (1, 3), (1, 4), (1, 5)]
    for bond in test_bonds:
        assert bond in bonds
    

def test_itp_top_atoms(bf4_itp_file: ItpFile):
    """
    Test for the _itp_top_name function.
    """
    atoms, bonds = _itp_top_atoms(bf4_itp_file)
    # now bonds are referred to atom index starting from 0
    test_bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    for bond in test_bonds:
        assert bond in bonds
    test_atoms = [
        ('B1', 'BF4', 1),
        ('F2', 'BF4', 1),
        ('F3', 'BF4', 1),
        ('F4', 'BF4', 1),
        ('F5', 'BF4', 1)
    ]
    for atom in test_atoms:
        assert atom in atoms


def test_itp_top(bf4_itp_fname: str, bf4_itp_file):
    name, atoms, bonds = itp_top(bf4_itp_fname)
    atoms_test, bonds_test = _itp_top_atoms(bf4_itp_file)
    assert _itp_top_name(bf4_itp_file) == name
    assert atoms_test == atoms
    assert bonds_test == bonds