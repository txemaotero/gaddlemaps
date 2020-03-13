"""
Test for the _components_top submodule.
"""

import os
import pytest

from gaddlemaps.components import MoleculeTop, AtomTop
from gaddlemaps.components._components_top import (_are_connected,
                                                   _find_connected_atoms)

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def bf4_itp_fname() -> str:
    """
    File name with the itp of a BF4 molecule.
    """
    return os.path.join(ACTUAL_PATH, '../../data/BF4_AA.itp')


@pytest.fixture
def bf4_molecule(bf4_itp_fname) -> MoleculeTop:
    """
    A BF4 MoleculeTop.
    """
    return MoleculeTop(bf4_itp_fname)


@pytest.fixture
def bmim_molecule() -> MoleculeTop:
    """
    A BF4 MoleculeTop.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.itp')
    return MoleculeTop(fname)

def test_single_atom_top():
    """
    Tests the initialization and attributes of an atom without bonded atoms.
    """
    atom = AtomTop('B1', 'BF4', 1, 0)
    assert atom.name == 'B1'
    assert atom.resname == 'BF4'
    assert atom.resid == 1
    assert atom.index == 0
    assert not atom.bonds
    assert hash(atom) == 0
    assert atom.residname == '1BF4'

    atom_compare = AtomTop('B1', 'BF4', 1, 0)
    assert atom == atom_compare
    assert atom is not atom_compare
    
    # Same index different resid
    atom_compare = AtomTop('B1', 'BF4', 5, 0)
    assert atom == atom_compare

    # Same resid different index
    atom_compare = AtomTop('B1', 'BF4', 1, 1)
    assert atom != atom_compare
    

def test_connect_atom_top():
    """
    Tests the connectivity related methods for AtomTop.

    This also test the _are_connected function and the _find_connected_atoms
    one.
    """
    # Atom0 connected with the rest
    atom0 = AtomTop('B1', 'BF4', 1, 0)
    atom1 = AtomTop('B2', 'BF4', 1, 1)
    atom2 = AtomTop('B3', 'BF4', 1, 2)
    atom3 = AtomTop('B4', 'BF4', 1, 3)

    atom0.connect(atom1)
    atom0.connect(atom2)
    atom0.connect(atom3)

    assert atom0.bonds == {1, 2, 3}
    assert atom0.closest_atoms() == [1, 2]
    assert atom0.closest_atoms(3) == [1, 2, 3]
    assert atom0.closest_atoms(0) == []

    assert atom1.bonds == {0}
    assert atom1.closest_atoms() == [0]
    assert atom2.bonds == {0}
    assert atom3.bonds == {0}

    connected_atoms = [atom0, atom1, atom2, atom3]
    connected_atoms_index = []
    _find_connected_atoms(connected_atoms, 0, connected_atoms_index)
    assert sorted(connected_atoms_index) == [0, 1, 2, 3]
    assert _are_connected(connected_atoms)

    # Not connected atom
    atom4 = AtomTop('B4', 'BF4', 1, 4)
    not_connected_atoms = [atom0, atom1, atom2, atom3, atom4]

    connected_atoms_index = []
    _find_connected_atoms(not_connected_atoms, 0, connected_atoms_index)
    assert sorted(connected_atoms_index) == [0, 1, 2, 3]

    connected_atoms_index = []
    _find_connected_atoms(not_connected_atoms, 4, connected_atoms_index)
    assert connected_atoms_index == [4]

    assert not _are_connected(not_connected_atoms)
    

class TestMoleculeTop:
    """
    Class to wrap all the MoleculeTop test.
    """
    def test_initialization(self, bf4_itp_fname: str):
        """
        Initialization.
        """
        with pytest.raises(ValueError, match=r".*.txt .*"):
            molecule = MoleculeTop('test.txt')
        with pytest.raises(ValueError, match=r".*.txt .*"):
            molecule = MoleculeTop('test.itp', file_format=".txt")

        molecule = MoleculeTop(bf4_itp_fname)
        assert molecule.ftop == bf4_itp_fname
        assert molecule.name == 'BF4'

        atoms_info = [
            ('B1', 'BF4', 1),
            ('F2', 'BF4', 1),
            ('F3', 'BF4', 1),
            ('F4', 'BF4', 1),
            ('F5', 'BF4', 1)
        ]
        atoms_test = [AtomTop(*info, i) for i, info in enumerate(atoms_info)]
        assert molecule.atoms == atoms_test

    def test_methods(self, bf4_molecule: MoleculeTop,
                     bmim_molecule: MoleculeTop):
        """
        Tests simple methods and magic-methods.
        """
        assert isinstance(bf4_molecule[0], AtomTop)
        assert bf4_molecule[0] == AtomTop('B1', 'BF4', 1, 0)
        assert len(bf4_molecule) == 5
        assert bf4_molecule == bf4_molecule
        assert bf4_molecule != 1
        assert bf4_molecule != bmim_molecule
        assert bf4_molecule.resname_len_list == [('BF4', 5)]
        assert bmim_molecule.resname_len_list == [('BMIM', 25)]
        for index, atom in enumerate(bf4_molecule):
            assert bf4_molecule.index(atom) == index
    
    def test_resnames(self, bf4_molecule: MoleculeTop):
        """
        Resnames access and setter.
        """
        assert bf4_molecule.resnames == ['BF4']
        with pytest.raises(ValueError, match='new_resnames.*'):
            bf4_molecule.resnames = 'test'  # type: ignore
        with pytest.raises(ValueError, match='Expected.*'):
            bf4_molecule.resnames = ['test', 'test2']

        bf4_molecule.resnames = ['BMIM']
        for atom in bf4_molecule:
            assert atom.resname == 'BMIM'
        assert bf4_molecule.resnames == ['BMIM']

        