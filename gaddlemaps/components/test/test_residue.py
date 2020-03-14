"""
Tests for the _residue submodule.
"""

import os
import pytest
from pathlib import Path

import numpy as np

from gaddlemaps import rotation_matrix
from gaddlemaps.parsers import GroFile
from gaddlemaps.components import Residue, AtomGro, MoleculeTop

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def atom_gro1() -> AtomGro:
    info = (1, 'BMIM', 'N1', 1, 4.668, 3.571, 8.232, -0.2489, 0.2514, 0.1046)
    return AtomGro(info)


@pytest.fixture
def atom_gro2() -> AtomGro:
    info = (1, 'BMIM', 'H2', 2, 5.668, 3.571, 8.232, -0.2489, 0.2514, 0.1046)
    return AtomGro(info)


@pytest.fixture
def atom_gro3() -> AtomGro:
    info = (1, 'BF4', 'H2', 2, 5.668, 3.571, 8.232)
    return AtomGro(info)


@pytest.fixture
def residue(atom_gro1: AtomGro, atom_gro2: AtomGro) -> Residue:
    return Residue([atom_gro1, atom_gro2])


@pytest.fixture
def bf4_mtop() -> MoleculeTop:
    """
    A BF4 MoleculeTop.
    """
    return MoleculeTop(os.path.join(ACTUAL_PATH, '../../data/BF4_AA.itp'))


class TestAtomGro:
    """
    Class to wrap AtomGro tests.
    """
    def test_initialization(self):
        """
        Initialization test.
        """
        atom = AtomGro((1, 'BMIM', 'N1', 1, 4.668, 3.571, 8.232, -0.2489,
                        0.2514, 0.1046))
        assert atom.resid == 1 
        assert atom.resname == 'BMIM'
        assert atom.name == 'N1'
        assert atom.atomid == 1
        assert np.isclose(atom.position, (4.668, 3.571, 8.232)).all()
        assert np.isclose(atom.velocity, (-0.2489, 0.2514, 0.1046)).all()

        atom = AtomGro((1, 'BMIM', 'N1', 1, 4.668, 3.571, 8.232))
        assert atom.velocity is None

    def test_add(self, atom_gro1: AtomGro, atom_gro2: AtomGro,
                 atom_gro3: AtomGro):
        """
        Test sum of atoms.
        """
        with pytest.raises(ValueError):
            residue = atom_gro1 + atom_gro3
        
        residue = atom_gro1 + atom_gro2
        assert isinstance(residue, Residue)
        assert len(residue) == 2
        assert residue[0] == atom_gro1
        assert residue[1] == atom_gro2
        assert residue.residname == atom_gro1.residname

    def test_methods(self, atom_gro1: AtomGro, atom_gro2: AtomGro,
                     atom_gro3: AtomGro):
        """
        Tests simple methods and magic-methods.
        """
        assert atom_gro1 == atom_gro1
        assert atom_gro1 != atom_gro2
        copy_1 = atom_gro1.copy()
        assert copy_1 == atom_gro1
        assert copy_1 is not atom_gro1

        assert atom_gro1.residname == '1BMIM'
        assert atom_gro2.residname == '1BMIM'
        assert atom_gro3.residname == '1BF4'

        assert atom_gro1.element == 'N'
        assert atom_gro2.element == 'H'
        assert atom_gro3.element == 'H'

    def test_atom_gro_line(self, atom_gro1: AtomGro, atom_gro3: AtomGro):
        """
        Tests gro_line method.
        """
        info = [1, 'BMIM', 'N1', 1, 4.668, 3.571, 8.232, -0.2489, 0.2514, 0.1046]
        assert atom_gro1.gro_line() == info
        line = '    1BMIM    N1    1   4.668   3.571   8.232 -0.2489  0.2514  0.1046'
        line_compare = atom_gro1.gro_line(parsed=False)
        assert isinstance(line_compare, str)
        assert line_compare.strip() == line.strip()
        assert len(atom_gro3.gro_line()) == 7


class TestResidue:
    """
    Wraps all the tests for Residue class.
    """
    def test_initialization(self, atom_gro1: AtomGro, atom_gro2: AtomGro,
                            atom_gro3: AtomGro):
        """
        Initialization test.
        """
        with pytest.raises(ValueError, match='Empty.*'):
            residue = Residue([])
        with pytest.raises(ValueError, match='The atoms.*'):
            residue = Residue([atom_gro1, atom_gro2, atom_gro3])

        residue = Residue([atom_gro1, atom_gro2])
        assert residue._atoms_gro == [atom_gro1, atom_gro2]

    def test_resid_name(self, residue: Residue):
        """
        Test resid and resname setter.
        """
        residue.resname = 'BF4'
        assert residue.resname == 'BF4'
        assert residue.residname == '1BF4'
        for atom in residue:
            assert atom.resname == 'BF4'
            assert atom.residname == '1BF4'

        residue.resid = 2
        assert residue.residname == '2BF4'
        assert residue.resid == 2
        for atom in residue:
            assert atom.resid == 2
            assert atom.residname == '2BF4'

    def test_methods(self, residue: Residue, atom_gro1: AtomGro):
        """
        Test some methods, magic-methods and attributes.
        """
        assert residue[0] == atom_gro1
        assert residue[0] is atom_gro1
        assert len(residue) == 2
        assert residue == residue
        residue_copy = residue.copy()
        assert residue == residue_copy
        assert residue is not residue_copy
        atoms = residue.atoms
        for at1, at2 in zip(residue, atoms):
            assert at1 == at2
            assert at1 is not at2

        assert residue.atoms_ids == [1 ,2]
        residue.atoms_ids = [4, 3]
        assert residue.atoms_ids == [4 ,3]
        with pytest.raises(TypeError):
            residue.atoms_ids = [1, 'hola']  # type: ignore
        with pytest.raises(IndexError):
            residue.atoms_ids = [1, 1, 2]

    def test_add(self, residue: Residue, atom_gro1: AtomGro,
                 atom_gro3: AtomGro):
        """
        Test the addition of Residue with other objects.
        """
        with pytest.raises(ValueError, match=r"To add an Ato.*"):
            new_residue = residue + atom_gro3

        residue_copy = residue.copy()
        residue_copy.resname = 'BF4'
        with pytest.raises(ValueError, match=r"To add two.*"):
            new_residue = residue + residue_copy

        with pytest.raises(TypeError):
            new_residue = residue + 1  # type: ignore

        new_residue = residue + atom_gro1
        assert len(new_residue) == 3
        assert new_residue.residname == atom_gro1.residname

        residue_copy = residue.copy()
        double_residue = residue + residue_copy
        assert len(double_residue) == 4
        assert double_residue.residname == residue.residname
        assert double_residue.residname == residue_copy.residname

    def test_position_velocities(self, residue: Residue):
        """
        Tests for atoms position and velocities getter and setter.
        """
        # Positions
        pos_test = np.array([
            [4.668, 3.571, 8.232],
            [5.668, 3.571, 8.232],
        ])
        assert np.isclose(pos_test, residue.atoms_positions).all()
        pos_change = np.random.random((2, 3))
        residue.atoms_positions = pos_change
        assert np.isclose(pos_change, residue.atoms_positions).all()
        for atom, pos in zip(residue, pos_change):
            assert np.isclose(pos, atom.position).all()

        # Geometric center
        center_test = np.mean(pos_change, axis=0)
        assert np.isclose(center_test, residue.geometric_center).all()
        assert center_test[0] == residue.x
        assert center_test[1] == residue.y
        assert center_test[2] == residue.z

        # Distance to zero
        assert ((center_test**2).sum())**.5 == residue.distance_to_zero

        # Distance to self
        assert residue.distance_to(residue) == 0
        # Distance to zero
        assert residue.distance_to([0, 0, 0]) == residue.distance_to_zero
        # Distance to zero with box vector
        box_vects = np.diag(residue.geometric_center)
        assert residue.distance_to([0, 0, 0],
                                   box_vects=box_vects) == 0
        # Distance to zero with box vector with inv
        inv_box = np.linalg.inv(box_vects)
        assert residue.distance_to([0, 0, 0], box_vects=inv_box,
                                   inv=True) == 0

        # Velocities
        vel_test = np.array([
            [-0.2489, 0.2514, 0.1046],
            [-0.2489, 0.2514, 0.1046]
        ])
        assert np.isclose(vel_test, residue.atoms_velocities).all()
        vel_change = np.random.random((2, 3))
        residue.atoms_velocities = vel_change
        assert np.isclose(vel_change, residue.atoms_velocities).all()
        for atom, vel in zip(residue, vel_change):
            assert np.isclose(vel, atom.velocity).all()

        residue.atoms_velocities = None
        assert residue.atoms_velocities is None
        for atom in residue:
            assert atom.velocity is None

    def test_rotate(self, residue: Residue):
        """
        Tests for simple rotations.
        """
        rot_mat = rotation_matrix([0, 0, 1], 2*np.pi)
        old_pos = residue.atoms_positions
        residue.rotate(rot_mat)
        assert np.isclose(old_pos, residue.atoms_positions).all()

        rot_mat = rotation_matrix([0, 0, 1], 0)
        residue.rotate(rot_mat)
        assert np.isclose(old_pos, residue.atoms_positions).all()

    def test_move(self, residue: Residue):
        """
        Tests for the moving methods.
        """
        # move_to
        move_to = [0, 1, 2]
        old_pos = residue.atoms_positions
        displ = move_to - residue.geometric_center
        residue.move_to(move_to)
        assert np.isclose(move_to, residue.geometric_center).all()
        assert np.isclose(old_pos + displ, residue.atoms_positions).all()

        # move
        move = [0, 1, 2]
        old_pos = residue.atoms_positions
        old_center = residue.geometric_center
        residue.move(move)
        assert np.isclose(move + old_center, residue.geometric_center).all()
        assert np.isclose(old_pos + move, residue.atoms_positions).all()

    def test_write_gro(self, residue: Residue, tmp_path: Path):
        """
        Test for writing gro files.
        """
        subdir = tmp_path / "residue_test"
        subdir.mkdir()
        fgrotmp = str(subdir / "test_write_gro.gro")
        residue.write_gro(fgrotmp)
        atoms = []
        with GroFile(fgrotmp) as fgro:
            assert fgro.natoms == 2
            for line in fgro:
                atoms.append(AtomGro(line))
        new_residue = Residue(atoms)
        assert new_residue == residue
        
    def test_update_from_top(self, residue: Residue, bf4_mtop: MoleculeTop):
        """
        Test for update_from_molecule_top method.
        """
        # Create a residue with 5 atoms
        new_residue = residue + residue + residue[0]
        assert len(new_residue) == 5
        names = ['N1', 'H2', 'N1', 'H2', 'N1']
        for atom, name in zip(new_residue, names):
            assert atom.name == name

        new_residue.update_from_molecule_top(bf4_mtop)
        names = ['B1', 'F2', 'F3', 'F4', 'F5'] 
        for atom, name in zip(new_residue, names):
            assert atom.name == name
