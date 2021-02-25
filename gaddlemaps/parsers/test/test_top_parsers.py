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
Tests for the _top_parsers submodule.
"""

import os
from pathlib import Path

import pytest

from gaddlemaps.parsers import ItpFile, read_topology
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


def test_itp_top(bf4_itp_fname: str, bf4_itp_file: ItpFile,
                 tmp_path: Path):
    name, atoms, bonds = read_topology(bf4_itp_fname)
    atoms_test, bonds_test = _itp_top_atoms(bf4_itp_file)
    assert _itp_top_name(bf4_itp_file) == name
    assert atoms_test == atoms
    assert bonds_test == bonds

    # Raises tests
    subdir = tmp_path / "itp_top_test"
    subdir.mkdir()
    fitptmp = str(subdir / "test_top_itp.gro")

    with open(fitptmp, 'w') as fitp:
        fitp.write('[ test ]')
    with pytest.raises(IOError, match=r'.*moleculetype.*'):
        _ = read_topology(fitptmp, file_format="itp")

    with open(fitptmp, 'w') as fitp:
        fitp.write('[ moleculetype ]\n')
    with pytest.raises(IOError, match=r'.*atoms.*'):
        _ = read_topology(fitptmp, file_format="itp")

    with open(fitptmp, 'w') as fitp:
        fitp.write('[ moleculetype ]\n[ atoms ]\n')
    with pytest.raises(IOError, match=r'There is not molecule.*'):
        _ = read_topology(fitptmp, file_format="itp")

    with open(fitptmp, 'w') as fitp:
        fitp.write('[ moleculetype ]\nTest   1\n[ atoms ]\n')
    with pytest.raises(IOError, match=r'There are not atoms.*'):
        _ = read_topology(fitptmp, file_format='itp')
