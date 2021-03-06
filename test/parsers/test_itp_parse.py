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
Tests for the _itp_parse submodule.
"""

import os
from pathlib import Path

import pytest

from gaddlemaps.parsers import (ItpFile, ItpLine, ItpLineAtom, ItpLineBonds,
                                ItpSection)

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

@pytest.fixture
def itp_file() -> ItpFile:
    data = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/SDS_AA.itp')
    return ItpFile(data)


@pytest.fixture
def dna_itp() -> ItpFile:
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_CG.itp')
    return ItpFile(fitp)


def test_itp_file(itp_file: ItpFile):
    """
    Test the ItpFile class.
    """
    keys = ('moleculetype', 'bonds', 'atoms', 'pairs')
    for key in keys:
        assert key in itp_file
        assert key in str(itp_file)

    assert itp_file is not itp_file.copy()
    assert 'itp file.' in repr(itp_file)


def test_itp_section(itp_file: ItpFile):
    """
    Test the ItpSection class.
    """
    header = itp_file['header']
    atoms = itp_file['atoms']
    bonds = itp_file['bonds']
    angles = itp_file['angles']
    assert bonds.section_name == 'bonds'
    # Instances
    assert isinstance(header, list)
    assert isinstance(atoms, list)
    assert isinstance(atoms, ItpSection)
    assert isinstance(angles, ItpSection)
    assert len(atoms) == 17
    assert isinstance(atoms[0], ItpLine)
    assert isinstance(atoms[0], ItpLineAtom)
    assert isinstance(bonds[0], ItpLineBonds)

    # Different initialization
    sec = ItpSection('test', ['line1', 'line2'])
    assert sec.section_name == 'test'
    assert isinstance(sec[0], ItpLine)
    assert repr(sec) == 'Itp "test" section.'
    for line, line_comp in zip(sec.lines, ['line1', 'line2']):
        assert str(line) == line_comp


def test_itp_line(itp_file: ItpFile):
    """
    Test the ItpFile class.
    """
    atom_line = itp_file['atoms'][0]
    # general line
    assert atom_line.content == '1	  S	1  SDS      S 	1     1.284	32.0600'
    assert atom_line.comment == ''
    assert atom_line.line == '     1	  S	1  SDS      S 	1     1.284	32.0600      \n'
    # Atom line
    assert atom_line.number == 1
    assert atom_line.type == 'S'
    assert atom_line.resid == 1
    assert atom_line.resname == 'SDS'
    assert atom_line.name == 'S'
    assert atom_line.cgnr == 1
    assert atom_line.charge == 1.284
    atom_line.charge = 1.5
    assert atom_line.charge == 1.5
    atom_line.charge = 1.284
    assert atom_line.mass == 32.06

    assert atom_line.parsed_line == [1, 'S', 1, 'SDS', 'S', 1, 1.284, 32.06]
    line = '1	  S	1  SDS      S 	1     1.284	32.0600'
    assert ItpLineAtom.read_itp_atom(line) == [1, 'S', 1, 'SDS', 'S', 1, 1.284,
                                               32.06]
    # Bond line
    bond_line = itp_file['bonds'][0]
    assert bond_line.atom_from == 17
    assert bond_line.atom_to == 16
    assert bond_line.funct == 1
    assert bond_line.const0 == 'gb_27'
    with pytest.raises(AttributeError):
        _ = bond_line.const1
    # Moleculetype line
    m_line = itp_file['moleculetype'][0]
    assert m_line.name == 'SDS'
    assert m_line.nrexcl == 3
    with pytest.raises(TypeError):
        m_line.name = 9
    with pytest.raises(ValueError):
        m_line.name = 'SDS Mod'
    with pytest.raises(ValueError):
        m_line.nrexcl = 0
    with pytest.raises(TypeError):
        m_line.nrexcl = 'hola'

    # Special cases
    test_str = '# test'
    line_itp = ItpLine(test_str)
    assert line_itp.line == test_str
    with pytest.warns(UserWarning):
        line_itp.content = 'content test'

    test_str = 'hello ; world'
    line_itp = ItpLine(test_str)
    assert 'hello' in line_itp.line
    assert ';' in line_itp.line
    assert 'world' in line_itp.line
    line_itp.comment = 'test comment'
    assert 'test comment' in line_itp.line

    with pytest.raises(IOError):
        _ = ItpLine('[ atoms ]')

    # no comment
    line_itp = ItpLine('test ;')
    assert line_itp.line.strip() == 'test'


def test_itp_line_atom():
    """
    Some test for atom lines
    """
    # Not full atom info
    line = '1	  S	1  SDS      S 	1'
    itp_line = ItpLineAtom(line)
    assert itp_line.charge is None
    with pytest.raises(AttributeError):
        _ = itp_line.mass
    itp_line.mass = 3
    assert itp_line.mass == 3

    with pytest.raises(AttributeError):
        _ = itp_line.pair_interaction
    itp_line.pair_interaction = 0
    assert itp_line.pair_interaction == 0

    with pytest.raises(AttributeError):
        _ = itp_line.aux_name
    itp_line.aux_name = 0
    assert itp_line.aux_name == 0


def test_itp_line_bond():
    """
    Some test for bond lines.
    """
    with pytest.raises(IOError):
        _ = ItpLineBonds('1')
    itp_line = ItpLineBonds('1  2')


def test_dna(dna_itp: ItpFile):
    """
    Test the ItpFile class for the DNA file.
    """
    assert 'moleculetype' in dna_itp
    assert 'atoms' in dna_itp
    assert 'bonds' in dna_itp
    assert 'constraints' in dna_itp
    assert 'angles' in dna_itp
    assert 'dihedrals' in dna_itp


def test_itp_write(itp_file: ItpFile, tmp_path: Path):
    subdir = tmp_path / "itp_test"
    subdir.mkdir()
    fitptmp = str(subdir / "test_compare_tmp.itp")
    itp_file.write(fitptmp)
    new_file = ItpFile(fitptmp)
    assert new_file.keys() == itp_file.keys()
    header = True
    for sec1, sec2 in zip(itp_file.values(), new_file.values()):
        assert len(sec1) == len(sec2)
        if header:
            print(type(sec1))
            assert isinstance(sec1, list)
            assert isinstance(sec2, list)
            header = False
        else:
            assert isinstance(sec1, ItpSection)
            assert isinstance(sec2, ItpSection)
            for l1, l2 in zip(sec1, sec2):
                print(l1, l2)
                assert l1.content == l2.content
                assert l1.comment == l2.comment
