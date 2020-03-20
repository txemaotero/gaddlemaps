'''
Tests for the GroFile class and functions defined in the init of the
submodule.
'''


import os
import pytest
import numpy as np
from gaddlemaps.parsers import (GroFile, _validate_res_atom_numbers,
                                extract_lattice_gro, dump_lattice_gro)


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
    assert open_fgro.natoms == 25
    box = open_fgro.box_matrix
    assert isinstance(box, np.ndarray)
    assert box.shape == (3, 3)
    box_test = np.array([
        [0.99, 0, 0],
        [0, 0.4, 0],
        [0, 0, 0.178]
    ])
    assert np.isclose(box, box_test).all()

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
