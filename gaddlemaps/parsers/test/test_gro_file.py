'''
Tests for the GroFile class.
'''


import os
import numpy as np
from gaddlemaps.parsers import GroFile


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


def test_str_gro_file():
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
    assert isinstance(first_line, list)
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
