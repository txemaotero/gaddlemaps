"""
Tests for the _itp_parse submodule.
"""

import os
import pytest
from pathlib import Path

from gaddlemaps.parsers import (ItpFile, ItpSection, ItpLineAtom, ItpLine,
                                ItpLineBonds)


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

@pytest.fixture
def itp_file() -> ItpFile:
    data = os.path.join(ACTUAL_PATH, '../../data/SDS_AA.itp')
    return ItpFile(data)


@pytest.fixture
def dna_itp() -> ItpFile:
    fitp = os.path.join(ACTUAL_PATH, '../../data/DNA_CG.itp')
    return ItpFile(fitp)


def test_itp_file(itp_file: ItpFile):
    """
    Test the ItpFile class.
    """
    assert 'moleculetype' in itp_file
    assert 'bonds' in itp_file
    assert 'atoms' in itp_file
    assert 'pairs' in itp_file
    assert itp_file is not itp_file.copy()


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
    assert atom_line.resnr == 1
    assert atom_line.resname == 'SDS'
    assert atom_line.atomname == 'S'
    assert atom_line.cgnr == 1
    assert atom_line.charge == 1.284
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
    
