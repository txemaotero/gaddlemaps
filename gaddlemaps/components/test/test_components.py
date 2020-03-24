"""
Tests for _components submodule.
"""

import os
from pathlib import Path
from typing import Callable, List

import numpy as np
import pytest

from gaddlemaps import rotation_matrix
from gaddlemaps.components import (Atom, AtomGro, AtomTop, Molecule,
                                   MoleculeTop, Residue, System)
from gaddlemaps.components._components import _molecule_top_and_residues_match
from gaddlemaps.parsers import GroFile

AtomGenerator = Callable[..., AtomGro]

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def atom_gro_generator() -> AtomGenerator:
    """
    A function to create AtomGro instances. 
    """
    def generate_atom_gro(resid: int = 1, resname: str = 'BMIM',
                          name: str = 'N1', atomid: int = 1) -> AtomGro:
        info = (resid, resname, name, atomid, 4.668,
                3.571, 8.232, -0.2489, 0.2514, 0.1046)
        return AtomGro(info)
    return generate_atom_gro

@pytest.fixture
def atom_n1_bmim(atom_gro_generator: AtomGenerator) -> Atom:
    """
    N1 Atom instance of BMIM.
    """
    atom_gro = atom_gro_generator()
    atom_top = AtomTop('N1', 'BMIM', 1, 0)
    return Atom(atom_top, atom_gro)


@pytest.fixture
def molecule_top_bmim() -> MoleculeTop:
    """
    BMIM MoleculeTop instance.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.itp')
    return MoleculeTop(fname)


@pytest.fixture
def residue_bmim() -> Residue:
    """
    BMIM Residue instance.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.gro')
    with GroFile(fname) as fgro:
        atoms = [AtomGro(line) for line in fgro]
    return Residue(atoms)


@pytest.fixture
def molecule_bmim(molecule_top_bmim: MoleculeTop,
                  residue_bmim: Residue) -> Molecule:
    """
    BMIM Molecule instance.
    """
    return Molecule(molecule_top_bmim, [residue_bmim])


@pytest.fixture
def molecule_top_protein() -> MoleculeTop:
    """
    MoleculeTop instance for molecule with multiple residues.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/Protein_CG.itp')
    return MoleculeTop(fname)


@pytest.fixture
def residue_protein() -> List[Residue]:
    """
    List of Residue instances for a molecule with multiple residues.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/Protein_CG.gro')
    residues: List[Residue] = []
    with GroFile(fname) as fgro:
        atom = AtomGro(next(fgro))
        residname = atom.residname
        atoms = [atom]
        for line in fgro:
            atom = AtomGro(line)
            if residname != atom.residname:
                residues.append(Residue(atoms))
                residname = atom.residname
                atoms = [atom]
            else:
                atoms.append(atom)
    residues.append(Residue(atoms))
    return residues


@pytest.fixture
def molecule_protein(molecule_top_protein: MoleculeTop,
                     residue_protein: List[Residue]) -> Molecule:
    """
    Molecule instance with multiple residues.
    """
    return Molecule(molecule_top_protein, residue_protein)


def test_molecule_top_and_residues_match(molecule_top_bmim: MoleculeTop,
                                         residue_bmim: Residue,
                                         molecule_top_protein: MoleculeTop,
                                         residue_protein: List[Residue]):
    """
    Test the match between residues and MoleculeTop.
    """
    # different length
    assert not _molecule_top_and_residues_match(molecule_top_bmim,
                                                residue_protein)
                                                
    assert _molecule_top_and_residues_match(molecule_top_bmim,
                                            [residue_bmim])
    assert _molecule_top_and_residues_match(molecule_top_protein,
                                            residue_protein)


class TestAtom:
    """
    Wraps all the tests for the Atom class.
    """
    def test_initialization(self, atom_gro_generator: AtomGenerator):
        """
        Tests Atom initialization and own attributes
        """
        atom_gro = atom_gro_generator()
        atom_top = AtomTop('N1', 'BMIM', 1, 0)
        atom = Atom(atom_top, atom_gro)
        assert atom._atom_gro is atom_gro
        assert atom._atom_top is atom_top
        new_atom = atom.copy()
        assert new_atom == atom
        assert new_atom is not atom
        assert hash(atom) == 0

        with pytest.raises(TypeError, match="atom_gro.*"):
            atom = Atom(atom_top, atom_top)  # type:ignore
        with pytest.raises(TypeError, match="atom_top.*"):
            atom = Atom(atom_gro, atom_gro)  # type:ignore

        atom_gro.name = 'B5'
        with pytest.raises(IOError):
            atom = Atom(atom_top, atom_gro)

        atom_gro.name = 'N1'
        atom_gro.resname = 'BF4'
        with pytest.raises(IOError):
            atom = Atom(atom_top, atom_gro)

    def test_get_set_attributes(self, atom_n1_bmim: Atom):
        """
        Tests get and set attributes for Atom.
        """
        # Name
        assert atom_n1_bmim.name == 'N1'
        atom_n1_bmim.name = 'B5'
        assert atom_n1_bmim.name == 'B5'
        assert atom_n1_bmim._atom_gro.name == 'B5'
        assert atom_n1_bmim._atom_top.name == 'B5'

        # Resname
        assert atom_n1_bmim.resname == 'BMIM'
        atom_n1_bmim.resname = 'BF4'
        assert atom_n1_bmim.resname == 'BF4'
        assert atom_n1_bmim._atom_gro.resname == 'BF4'
        assert atom_n1_bmim._atom_top.resname == 'BF4'

        # resids
        with pytest.raises(AttributeError):
            resid = atom_n1_bmim.resid
        with pytest.raises(AttributeError):
            atom_n1_bmim.resid = 1

        assert atom_n1_bmim.gro_resid == 1
        atom_n1_bmim.gro_resid = 3
        assert atom_n1_bmim.gro_resid == 3
        assert atom_n1_bmim._atom_gro.resid == 3

        assert atom_n1_bmim.top_resid == 1
        atom_n1_bmim.top_resid = 5
        assert atom_n1_bmim.top_resid == 5
        assert atom_n1_bmim._atom_top.resid == 5

        # Residname is the gro one
        assert atom_n1_bmim.residname == '3BF4'

        # index
        assert atom_n1_bmim.index == 0
        assert hash(atom_n1_bmim) == 0
        atom_n1_bmim.index = 5
        assert atom_n1_bmim.index == 5
        assert hash(atom_n1_bmim) == 5
        assert atom_n1_bmim._atom_top.index == 5
        
        # bonds 
        assert atom_n1_bmim.bonds == atom_n1_bmim._atom_top.bonds
        atom_n1_bmim.bonds = {1, 32, 4}
        assert atom_n1_bmim.bonds == {1, 32, 4}
        assert atom_n1_bmim.closest_atoms() == [1, 4]
        atom_n1_bmim.connect(atom_n1_bmim)
        assert atom_n1_bmim.bonds == {1, 32, 4, 5}

        # atomid
        assert atom_n1_bmim.atomid == 1
        atom_n1_bmim.atomid = 6
        assert atom_n1_bmim.atomid == 6
        assert atom_n1_bmim._atom_gro.atomid == 6

        # element
        assert atom_n1_bmim.element == 'B'

        # position
        assert np.isclose(atom_n1_bmim.position, [4.668, 3.571, 8.232]).all()
        atom_n1_bmim.position = np.array([0, 0, 0])
        assert np.isclose(atom_n1_bmim.position, [0, 0, 0]).all()
        assert np.isclose(atom_n1_bmim._atom_gro.position, [0, 0, 0]).all()

        # velocity
        assert np.isclose(atom_n1_bmim.velocity, [-0.2489, 0.2514, 0.1046]).all()
        atom_n1_bmim.velocity = None
        assert atom_n1_bmim.velocity is None
        assert atom_n1_bmim._atom_gro.velocity is None

        # gro_line
        assert atom_n1_bmim.gro_line() == atom_n1_bmim._atom_gro.gro_line()


class TestMolecule:
    """
    Wraps all the tests for the Molecule class.
    """
    def test_init_one_residue(self, molecule_top_bmim: MoleculeTop,
                              residue_bmim: Residue):
        """
        Tests Molecule initialization with only one residue and some basic
        attributes and methods.
        """
        molec = Molecule(molecule_top_bmim, [residue_bmim])

        assert len(molec) == len(molecule_top_bmim)
        assert len(molec) == len(residue_bmim)
        assert len(molec) == 25

        for atom, atom_top, atom_gro in zip(molec, molecule_top_bmim,
                                            residue_bmim):
            atom_compare = Atom(atom_top, atom_gro)
            assert atom == atom_compare

        for atom, atom_comp in zip(molec, molec.atoms):
            assert atom == atom_comp

        assert molec[-1] == atom_compare
        assert atom_compare != molec[0]
        assert molec.index(atom_compare) == 24

        assert molec._each_atom_resid == [0] * len(molec)

        assert molec == molec
        assert molec is molec

        assert molec.copy() == molec
        assert molec.copy() is not molec

        assert len(molec.residues) == 1
        assert molec.residues[0] == residue_bmim
        assert molec.residues[0] is not residue_bmim
        assert molec.molecule_top is molecule_top_bmim

    def test_init_multi_residue(self, molecule_top_protein: MoleculeTop,
                                residue_protein: List[Residue]):
        """
        Tests Molecule initialization with multiple residues and some basic
        attributes and methods.
        """
        molec = Molecule(molecule_top_protein, residue_protein)

        assert len(molec) == len(molecule_top_protein)
        assert len(molec) == sum(len(res) for res in residue_protein)
        assert len(molec) == 34

        index = 0
        for res in residue_protein:
            for atom_gro in res:
                atom_compare = Atom(molecule_top_protein[index], atom_gro)
                assert molec[index] == atom_compare
                index += 1

        for atom, atom_comp in zip(molec, molec.atoms):
            assert atom == atom_comp

        assert molec[-1] == atom_compare
        assert molec.index(atom_compare) == 33

        assert molec._each_atom_resid == [0, 0, 0, 1, 1, 1, 2, 2, 3,
                                          3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7,
                                          8, 8, 8, 9, 9, 10, 10, 10, 11, 11,
                                          11, 12, 12]

        assert molec == molec
        assert molec is molec

        assert molec.copy() == molec
        assert molec.copy() is not molec

        assert len(molec.residues) == 13
        for res1, res2 in zip(molec.residues, residue_protein):
            assert res1 == res2
            assert res1 is not res2
        assert molec.molecule_top is molecule_top_protein

    def test_init_from_files(self):
        """
        Test the molecule initialization from files.
        """
        fitp = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.itp')
        fgro = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.gro')
        bmim = Molecule.from_files(fgro, fitp)
        assert len(bmim) == 25
        fitp = os.path.join(ACTUAL_PATH, '../../data/BMIM_CG.itp')
        fgro = os.path.join(ACTUAL_PATH, '../../data/system_bmimbf4_cg.gro')
        with pytest.raises(IOError):
            bmim = Molecule.from_files(fgro, fitp)

    def test_resnames(self, molecule_protein: Molecule,
                      molecule_bmim: Molecule):
        """
        Tests for the resnames property
        """
        assert molecule_bmim.resnames == ['BMIM']
        molecule_bmim.resnames = ['BF4']
        assert molecule_bmim.resnames == ['BF4']
        assert molecule_bmim.molecule_top.resnames == ['BF4']
        assert molecule_bmim.residues[0].resname == 'BF4'
        for atom in molecule_bmim:
            assert atom.resname == 'BF4'

        molecule_bmim.resnames = 'test'  # type: ignore
        assert molecule_bmim.resnames == ['test']
        assert molecule_bmim.molecule_top.resnames == ['test']
        assert molecule_bmim.residues[0].resname == 'test'
        for atom in molecule_bmim:
            assert atom.resname == 'test'

        resnames_compare = ['ARG', 'ARG', 'LEU', 'LYS', 'ARG', 'LEU', 'LEU',
                            'ARG', 'ARG', 'LEU', 'LYS', 'ARG', 'LEU']
        assert molecule_protein.resnames == resnames_compare
        resnames_compare[3] = 'Test'
        molecule_protein.resnames = resnames_compare
        assert molecule_protein.resnames == resnames_compare
        for res, res_test in zip(molecule_protein.residues, resnames_compare):
            assert res.resname == res_test
        for index, atom in enumerate(molecule_protein):
            res_index = molecule_protein._each_atom_resid[index]
            assert atom.resname == resnames_compare[res_index]

        molecule_protein.resnames = 'hello'  # type: ignore
        assert molecule_protein.resnames == ['hello'] * len(resnames_compare)
        for res, res_test in zip(molecule_protein.residues, resnames_compare):
            assert res.resname == 'hello'
        for atom in molecule_protein:
            assert atom.resname == 'hello'

    def test_resids(self, molecule_protein: Molecule,
                    molecule_bmim: Molecule):
        """
        Tests for the resids property
        """
        assert molecule_bmim.resids == [1]
        molecule_bmim.resids = [5]
        assert molecule_bmim.resids == [5]
        assert molecule_bmim.molecule_top.resids == [5]
        assert molecule_bmim.residues[0].resid == 5
        for atom in molecule_bmim:
            assert atom.top_resid == 5
            assert atom.gro_resid == 5

        molecule_bmim.resids = 3  # type: ignore
        assert molecule_bmim.resids == [3]
        assert molecule_bmim.molecule_top.resids == [3]
        assert molecule_bmim.residues[0].resid == 3
        for atom in molecule_bmim:
            assert atom.top_resid == 3
            assert atom.gro_resid == 3

        resids_compare = list(range(1, 14))
        assert molecule_protein.resids == resids_compare
        resids_compare[3] = 8
        molecule_protein.resids = resids_compare
        assert molecule_protein.resids == resids_compare
        for res, res_test in zip(molecule_protein.residues, resids_compare):
            assert res.resid == res_test
        for index, atom in enumerate(molecule_protein):
            res_index = molecule_protein._each_atom_resid[index]
            assert atom.top_resid == resids_compare[res_index]
            assert atom.gro_resid == resids_compare[res_index]

        molecule_protein.resids = 49  # type: ignore
        assert molecule_protein.resids == [49] * len(resids_compare)
        for res, res_test in zip(molecule_protein.residues, resids_compare):
            assert res.resid == 49
        for atom in molecule_protein:
            assert atom.top_resid == 49
            assert atom.gro_resid == 49

    def test_get_set_attributes(self, molecule_top_bmim: MoleculeTop,
                                molecule_bmim: Molecule):
        """
        Tests get and set attributes for Molecule. Molecule getattr just access
        the molecule_top attributes.
        """
        # Name
        assert molecule_top_bmim.name == molecule_bmim.name
        molecule_bmim.name = 'Test'
        assert molecule_bmim.name == 'Test'
        assert molecule_bmim.molecule_top.name == 'Test'

        # Resname len list
        assert molecule_bmim.resname_len_list == [('BMIM', 25)]

    def test_bonds_distance(self, molecule_bmim: Molecule):
        """
        Test bonds distances for bmim.
        """
        bonds_dist = molecule_bmim.bonds_distance
        for index, atom in enumerate(molecule_bmim):
            assert set(val[0] for val in bonds_dist[index]) == atom.bonds
            for index_to, dist in bonds_dist[index]:
                diff = atom.position - molecule_bmim[index_to].position
                dist_test = np.linalg.norm(diff)
                assert round(dist, 5) == round(dist_test, 5)

    def test_position_velocities_bmim(self, molecule_bmim: Molecule):
        """
        Test for the atom position and velocities get and setter for bmim
        molecule.
        """
        test_bmim_pos = np.array([
            [1.593, 1.896, 0.729],
            [1.706, 1.984, 0.708],
            [1.807, 1.892, 0.660],
            [1.755, 1.761, 0.652],
            [1.626, 1.764, 0.693],
            [1.738, 2.029, 0.805],
            [1.809, 1.674, 0.618],
            [1.559, 1.679, 0.698],
            [1.464, 1.937, 0.778],
            [1.942, 1.930, 0.624],
            [2.044, 1.858, 0.713],
            [2.188, 1.899, 0.677],
            [2.289, 1.828, 0.765],
            [1.463, 2.046, 0.799],
            [1.441, 1.883, 0.872],
            [1.385, 1.914, 0.703],
            [1.960, 1.904, 0.517],
            [1.955, 2.039, 0.636],
            [2.024, 1.884, 0.819],
            [2.034, 1.748, 0.701],
            [2.200, 2.009, 0.689],
            [2.208, 1.872, 0.571],
            [2.392, 1.859, 0.737],
            [2.273, 1.854, 0.872],
            [2.281, 1.718, 0.753]
        ])
        assert np.isclose(test_bmim_pos, molecule_bmim.atoms_positions).all()
        assert molecule_bmim.atoms_velocities is None

        test_bmim_pos[1] = np.array((1, 1, 1))
        molecule_bmim.atoms_positions = test_bmim_pos
        assert np.isclose(test_bmim_pos, molecule_bmim.atoms_positions).all()

        change_atom = np.array((2, 3, 4))
        test_bmim_pos[3] = change_atom
        molecule_bmim[3].position = change_atom
        assert np.isclose(molecule_bmim[3].position, change_atom).all()
        assert np.isclose(test_bmim_pos, molecule_bmim.atoms_positions).all()

        geom_test = np.mean(test_bmim_pos, axis=0)
        assert np.isclose(geom_test, molecule_bmim.geometric_center).all()

        with pytest.raises(ValueError):
            molecule_bmim.atoms_positions = change_atom

        with pytest.raises(ValueError):
            molecule_bmim.atoms_velocities = change_atom

        molecule_bmim.atoms_velocities = test_bmim_pos
        assert np.isclose(test_bmim_pos, molecule_bmim.atoms_velocities).all()

        change_atom = np.array((5, 1, 3))
        test_bmim_pos[5] = change_atom
        assert np.isclose(molecule_bmim[5].velocity, change_atom).all()
        assert np.isclose(test_bmim_pos, molecule_bmim.atoms_velocities).all()

    def test_position_velocities_protein(self, molecule_protein: Molecule):
        """
        Test for the atom position and velocities get and setter for the
        molecule with multiple residues.
        """
        test_protein_pos = np.array([
            [2.490, 3.850, 5.330],
            [2.310, 4.040, 5.170],
            [2.190, 4.380, 5.220],
            [2.480, 3.530, 5.200],
            [2.280, 3.350, 5.350],
            [2.180, 3.040, 5.290],
            [2.700, 3.410, 5.120],
            [2.770, 3.370, 4.800],
            [2.900, 3.290, 5.390],
            [3.220, 3.190, 5.390],
            [3.440, 3.130, 5.540],
            [2.700, 3.030, 5.470],
            [2.490, 3.060, 5.730],
            [2.430, 3.160, 6.050],
            [2.710, 2.860, 5.200],
            [2.550, 2.950, 4.940],
            [3.040, 2.820, 5.150],
            [3.110, 2.980, 4.860],
            [3.090, 2.690, 5.460],
            [3.280, 2.580, 5.720],
            [3.120, 2.490, 5.990],
            [2.830, 2.450, 5.450],
            [2.750, 2.270, 5.670],
            [2.460, 2.210, 5.800],
            [2.910, 2.280, 5.170],
            [2.850, 2.430, 4.850],
            [3.210, 2.310, 5.160],
            [3.250, 2.010, 5.010],
            [3.330, 1.820, 5.160],
            [3.220, 2.140, 5.500],
            [3.550, 2.070, 5.510],
            [3.830, 2.170, 5.540],
            [3.050, 1.870, 5.590],
            [2.790, 1.770, 5.390]
        ])
        assert np.isclose(test_protein_pos, molecule_protein.atoms_positions).all()
        assert molecule_protein.atoms_velocities is None

        test_protein_pos[1] = np.array((1, 1, 1))
        molecule_protein.atoms_positions = test_protein_pos
        assert np.isclose(test_protein_pos, molecule_protein.atoms_positions).all()

        change_atom = np.array((2, 3, 4))
        test_protein_pos[3] = change_atom
        molecule_protein[3].position = change_atom
        assert np.isclose(molecule_protein[3].position, change_atom).all()
        assert np.isclose(test_protein_pos, molecule_protein.atoms_positions).all()

        geom_test = np.mean(test_protein_pos, axis=0)
        assert np.isclose(geom_test, molecule_protein.geometric_center).all()

        with pytest.raises(ValueError):
            molecule_protein.atoms_positions = change_atom

        with pytest.raises(ValueError):
            molecule_protein.atoms_velocities = change_atom

        molecule_protein.atoms_velocities = test_protein_pos
        assert np.isclose(test_protein_pos, molecule_protein.atoms_velocities).all()

        change_atom = np.array((5, 1, 3))
        test_protein_pos[5] = change_atom
        assert np.isclose(molecule_protein[5].velocity, change_atom).all()
        assert np.isclose(test_protein_pos, molecule_protein.atoms_velocities).all()

    def test_atoms_ids(self, molecule_protein: Molecule):
        """
        Test for the atoms_ids property for the molecule with multiple residues.
        """
        test_ids = list(range(1, 35))
        assert molecule_protein.atoms_ids == test_ids

        test_ids[2] = 5
        molecule_protein.atoms_ids = test_ids
        assert molecule_protein.atoms_ids == test_ids

        test_ids[4] = 40
        molecule_protein[4].atomid = 40
        assert molecule_protein.atoms_ids == test_ids

        with pytest.raises(IndexError):
            molecule_protein.atoms_ids = [1, 2]

        with pytest.raises(TypeError):
            test_wrong = [1] * (len(test_ids) - 1) + ['test']  # type: ignore
            molecule_protein.atoms_ids = test_wrong

    def test_hidden_methods(self, molecule_bmim: Molecule):
        """
        Test not implemented properties and methods.
        """
        with pytest.raises(AttributeError):
            test = molecule_bmim.resname
        with pytest.raises(AttributeError):
            molecule_bmim.resname = 'test'

        with pytest.raises(AttributeError):
            test = molecule_bmim.resid # type: ignore
        with pytest.raises(AttributeError):
            molecule_bmim.resid = 'test'  # type: ignore

        with pytest.raises(AttributeError):
            test = molecule_bmim.residname
        with pytest.raises(AttributeError):
            molecule_bmim.residname = 'test'  # type: ignore

        with pytest.raises(AttributeError):
            test = molecule_bmim.remove_atom(molecule_bmim[0]) # type: ignore
        
    def test_rotate(self, molecule_protein: Molecule):
        """
        Test for the rotate method.
        """
        rot_mat = rotation_matrix([0, 0, 1], 2*np.pi)
        old_pos = molecule_protein.atoms_positions
        molecule_protein.rotate(rot_mat)
        assert np.isclose(old_pos, molecule_protein.atoms_positions).all()

        rot_mat = rotation_matrix([0, 0, 1], 0)
        molecule_protein.rotate(rot_mat)
        assert np.isclose(old_pos, molecule_protein.atoms_positions).all()
        
    def test_move(self, molecule_protein: Molecule):
        """
        Tests for the moving methods.
        """
        # move_to
        move_to = [0, 1, 2]
        old_pos = molecule_protein.atoms_positions
        displ = move_to - molecule_protein.geometric_center
        molecule_protein.move_to(move_to)
        assert np.isclose(move_to, molecule_protein.geometric_center).all()
        assert np.isclose(old_pos + displ, molecule_protein.atoms_positions).all()

        # move
        move = [0, 1, 2]
        old_pos = molecule_protein.atoms_positions
        old_center = molecule_protein.geometric_center
        molecule_protein.move(move)
        assert np.isclose(move + old_center, molecule_protein.geometric_center).all()
        assert np.isclose(old_pos + move, molecule_protein.atoms_positions).all()

    def test_write_gro(self, molecule_protein: Molecule, tmp_path: Path):
        """
        Test for writing gro files.
        """
        subdir = tmp_path / "molecule_test"
        subdir.mkdir()
        fgrotmp = str(subdir / "test_write_gro.gro")
        molecule_protein.write_gro(fgrotmp)
        with GroFile(fgrotmp) as fgro:
            assert fgro.natoms == 34
            for line, atom in zip(fgro, molecule_protein):
                atom_gro = AtomGro(line)
                assert atom.atom_gro == atom_gro

    def test_update_from_top(self, molecule_bmim: Molecule):
        """
        Test for update_from_molecule_top method.
        """
        # Make a copy of the original top
        fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.itp')
        mol_top = MoleculeTop(fname)
        for at1, at2 in zip(mol_top, molecule_bmim):
            assert at1.name == at2.name

        mol_top[1].name = 'test'
        molecule_bmim.update_from_molecule_top(mol_top)
        for at1, at2 in zip(mol_top, molecule_bmim):
            assert at1.name == at2.name
        assert molecule_bmim[1].name == 'test'

    def test_distance_to(self, molecule_protein: Molecule):
        """
        Test for method related with distances.
        """
        center = molecule_protein.geometric_center
        # Distance to zero
        assert ((center**2).sum())**.5 == molecule_protein.distance_to_zero

        # Distance to self
        assert molecule_protein.distance_to(molecule_protein) == 0
        # Distance to zero
        assert molecule_protein.distance_to([0, 0, 0]) == molecule_protein.distance_to_zero
        # Distance to zero with box vector
        box_vects = np.diag(molecule_protein.geometric_center)
        assert round(molecule_protein.distance_to([0, 0, 0],
                                                  box_vects=box_vects), 8) == 0
        # Distance to zero with box vector with inv
        inv_box = np.linalg.inv(box_vects)
        assert round(molecule_protein.distance_to([0, 0, 0], box_vects=inv_box,
                                                  inv=True), 8) == 0

    def test_copy(self, molecule_bmim: Molecule):
        """
        Test the copy method of Molecule class.
        """
        copy_mol = molecule_bmim.copy()

        assert copy_mol == molecule_bmim
        assert copy_mol.resnames == molecule_bmim.resnames

        assert copy_mol.molecule_top == molecule_bmim.molecule_top
        assert copy_mol.molecule_top is molecule_bmim.molecule_top

        assert copy_mol.residues == molecule_bmim.residues
        assert copy_mol.residues is not molecule_bmim.residues

        # Tricky change in the resname of the topology
        copy_mol.molecule_top.resnames = ['Test']
        assert molecule_bmim.molecule_top.resnames == ['Test']

    def test_deep_copy(self, molecule_bmim: Molecule):
        """
        Test the deep_copy method of Molecule class.
        """
        copy_mol = molecule_bmim.deep_copy()

        assert copy_mol == molecule_bmim
        assert copy_mol.resnames == molecule_bmim.resnames

        assert copy_mol.molecule_top == molecule_bmim.molecule_top
        assert copy_mol.molecule_top is not molecule_bmim.molecule_top

        assert copy_mol.residues == molecule_bmim.residues
        assert copy_mol.residues is not molecule_bmim.residues

        # Tricky change in the resname of the topology
        copy_mol.molecule_top.resnames = ['Test']
        assert molecule_bmim.molecule_top.resnames == ['BMIM']
