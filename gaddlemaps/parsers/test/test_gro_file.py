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
Tests for the GroFile class and functions defined in the init of the
submodule.
'''


import os
from pathlib import Path

import numpy as np
import pytest

from gaddlemaps.parsers import (GroFile, _validate_res_atom_numbers,
                                dump_lattice_gro, extract_lattice_gro)

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


def test_validate_res_atom_numbers():
    """
    Test for _validate_res_atom_numbers function.
    """
    line_test = '    1BMIM    C2    2   1.706   1.984   0.708'
    assert _validate_res_atom_numbers(line_test) == (1, 2)
    line_test = '    ABMIM    C2    2   1.706   1.984   0.708'
    with pytest.raises(IOError, match=r".* residue .*"):
        _validate_res_atom_numbers(line_test)
    line_test = '    1BMIM    C2    A   1.706   1.984   0.708'
    with pytest.raises(IOError, match=r".* atom .*"):
        _validate_res_atom_numbers(line_test)


def test_lattice_gro():
    line_test = '   0.99000   0.40000   0.17800'
    box = extract_lattice_gro(line_test)
    box_test = np.array([
        [0.99, 0, 0],
        [0, 0.4, 0],
        [0, 0, 0.178]
    ])
    assert np.isclose(box, box_test).all()
    assert line_test.strip() == dump_lattice_gro(box).strip()


def test_gro_file():
    """
    Test for opening a gro from string
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.gro')
    open_fgro = GroFile(fname, 'r')
    assert open_fgro.name == fname

    assert open_fgro.natoms == 25
    with pytest.raises(AttributeError):
        open_fgro.natoms = 4

    with pytest.raises(AttributeError):
        open_fgro.comment = 'Test line'

    with pytest.raises(IndexError):
        open_fgro.seek_atom(40)

    assert isinstance(open_fgro.position_format, tuple)
    assert len(open_fgro.position_format) == 2
    assert isinstance(open_fgro.position_format[0], int)
    assert isinstance(open_fgro.position_format[1], int)

    with pytest.raises(AttributeError):
        open_fgro.position_format = (1, 3)

    with pytest.raises(AttributeError):
        open_fgro.box_matrix = np.eye(3)

    box = open_fgro.box_matrix
    assert isinstance(box, np.ndarray)
    assert box.shape == (3, 3)
    box_test = np.array([
        [0.99, 0, 0],
        [0, 0.4, 0],
        [0, 0, 0.178]
    ])
    assert np.isclose(box, box_test).all()

    open_fgro.seek_atom(0)
    first_line = next(open_fgro)
    assert isinstance(first_line, tuple)
    assert len(first_line) == 7
    assert first_line[0] == 1
    assert first_line[1] == 'BMIM'
    assert first_line[2] == 'N1'
    assert first_line[3] == 1
    assert first_line[4] == 1.593
    assert first_line[5] == 1.896
    assert first_line[6] == 0.729
    # Second line not parsed
    second_line = open_fgro.readline(parsed=False)
    assert isinstance(second_line, str)
    test_second = '1BMIM    C2    2   1.706   1.984   0.708'
    assert second_line.strip() == test_second


def test_empty_gro_file(tmp_path: Path):
    """
    Test for opening an empty gro and some read write operations.
    """
    subdir = tmp_path / "gro_file_test"
    subdir.mkdir()
    fgrotmp = str(subdir / "test_write_fgro.gro")
    open_fgro = GroFile(fgrotmp, 'w')
    with pytest.raises(ValueError):
        _ = open_fgro.natoms
    with pytest.raises(ValueError):
        open_fgro.seek_atom(3)
    open_fgro.natoms = 1

    assert isinstance(open_fgro.position_format, tuple)
    assert len(open_fgro.position_format) == 2
    assert isinstance(open_fgro.position_format[0], int)
    assert isinstance(open_fgro.position_format[1], int)

    orig_format = open_fgro.position_format
    open_fgro.position_format = (1, 3)
    open_fgro.position_format = orig_format

    open_fgro.box_matrix = np.array((1, 1, 1))
    assert np.isclose(open_fgro.box_matrix, np.eye(3)).all()
    with pytest.raises(ValueError):
        open_fgro.box_matrix = np.array((1, 1, 1, 1))

    bad_fgro = str(subdir / "test_bad_fgro.gro")
    with open(bad_fgro, 'w') as fgro:
        fgro.write('Comment\nbad\n')
    with pytest.raises(IOError, match=r'Second.*'):
        _ = GroFile(bad_fgro)

    with open(bad_fgro, 'w') as fgro:
        fgro.write('Comment\n1\n')
        fgro.write('    1BMIM    C2    A   1.706   1.984   0.708\n')
        fgro.write('    1BMIM    C2    A   1.706   1.984   0.708\n')
        line_test = '   0.99000   0.40000   0.17800'
    with pytest.raises(IOError, match=r'The number of atoms.*'):
        _ = GroFile(bad_fgro)

    with open(bad_fgro, 'w') as fgro:
        fgro.write('Comment\n3\n')
        fgro.write('    1BMIM    C2    A   1.706   1.984   0.708\n')
        fgro.write('    1BMIM    C2    A   1.706   1.984   0.708\n')
        line_test = '   0.99000   0.40000   0.17800'
    with pytest.raises(IOError, match=r'The number of atoms.*'):
        _ = GroFile(bad_fgro)

    del open_fgro

    open_fgro = GroFile(fgrotmp, 'w')
    line = '    1BMIM    C2    1   1.706   1.984   0.708'
    open_fgro.writelines([line, line])  # type: ignore
    open_fgro.close()

    open_fgro = GroFile(fgrotmp)
    assert open_fgro.natoms == 2
    lines = open_fgro.readlines()
    assert isinstance(lines, list)
    assert isinstance(lines[0], tuple)
    open_fgro.close()

    test_list = [1, 2]
    with pytest.raises(ValueError):
        _ = GroFile.parse_atomlist(test_list)  # type: ignore
    with pytest.raises(IOError):
        _ = GroFile.parse_atomline('Test')

    with pytest.warns(UserWarning):
        open_fgro = GroFile(fgrotmp, 'w')
        open_fgro.close()

    open_fgro = GroFile(fgrotmp, 'w')
    open_fgro.natoms = 3
    line = '    1BMIM    C2    1   1.706   1.984   0.708'
    open_fgro.writelines([line, line])  # type: ignore
    with pytest.raises(IOError):
        open_fgro.close()

    with pytest.raises(ValueError):
        GroFile.determine_format('\n\n')
    with pytest.raises(IOError, match=r"Found.*"):
        GroFile.determine_format(line + '.\n')

    line = '    1BMIM    C2    1   17.00980986   1.984   0.7080'
    with pytest.raises(IOError, match=r"Some.*"):
        GroFile.determine_format(line + '\n')

    with pytest.warns(UserWarning):
        new = GroFile.validate_string('TestLong')
        assert new == 'TestL'