"""
Tests for the _system submodule.
"""

import os
from typing import Tuple

import numpy as np
import pytest

from gaddlemaps.components import MoleculeTop, System, SystemGro

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def bmimbf4_gro_fname() -> str:
    return os.path.join(ACTUAL_PATH, '../../data/system_bmimbf4_cg.gro')


@pytest.fixture
def system() -> System:
    fgro = os.path.join(ACTUAL_PATH, '../../data/system_CG.gro')
    fitpDNA = os.path.join(ACTUAL_PATH, '../../data/DNA_CG.itp')
    fitpDPSM = os.path.join(ACTUAL_PATH, '../../data/DPSM_CG.itp')
    fitpVTE = os.path.join(ACTUAL_PATH, '../../data/vitamin_E_CG.itp')
    return System(fgro, fitpDNA, fitpDPSM, fitpVTE)


@pytest.fixture
def system_bmimbf4_gro(bmimbf4_gro_fname: str) -> SystemGro:
    return SystemGro(bmimbf4_gro_fname)


class TestSystemGro:
    """
    Wraps all the tests for the SystemGro object.
    """
    def test_initialization(self, bmimbf4_gro_fname: str):
        """
        Test systemGro initialization with bmim bf4 system.
        """
        system = SystemGro(bmimbf4_gro_fname)
        assert system.fgro == bmimbf4_gro_fname
        assert len(system.different_molecules) == 2

        bf4 = system.different_molecules[0]
        assert bf4.resname == 'BF4'
        assert len(bf4) == 1

        bmim = system.different_molecules[1]
        assert bmim.resname == 'BMIM'
        assert len(bmim) == 3

        assert len(system) == 600
        atom_index = 1
        for index, res in enumerate(system):
            assert res.resid == index + 1
            if index < 300:
                assert res.resname == 'BF4'
            else:
                assert res.resname == 'BMIM'
            for atom in res:
                assert atom.atomid == atom_index
                atom_index += 1
        
        # Check the last positions and velocities
        test_positions = np.array([
            [0.288, 1.338, 0.925],
            [0.352, 1.069, 0.807],
            [0.540, 0.960, 0.772]
        ])
        assert np.isclose(res.atoms_positions, test_positions).all()
        test_velocities = np.array([
            [-0.1079, 0.1208, 0.0298],
            [0.1199, -0.0489, 0.5257],
            [0.0236, 0.0129, -0.2678]
        ])
        assert np.isclose(res.atoms_velocities, test_velocities).all()

        # Random mols positions and velocities
        test_positions = np.array([[2.463, 0.308, 3.294]])
        assert np.isclose(system[4].atoms_positions, test_positions).all()
        test_velocities = np.array([[0.0639, 0.1574, -0.0145]])
        assert np.isclose(system[9].atoms_velocities, test_velocities).all()

        test_positions = np.array([
            [4.031, 3.087, 1.147],
            [4.317, 3.012, 1.097],
            [4.098, 3.018, 1.102]
        ])
        assert np.isclose(system[420].atoms_positions, test_positions).all()

        test_velocities = np.array([
            [0.1281, -0.0416, 0.0584],
            [0.2010, -0.1297, -0.1044],
            [0.1105, 0.3189, 0.0302]
        ])
        assert np.isclose(system[366].atoms_velocities, test_velocities).all()

    def test_attributes(self, system_bmimbf4_gro: SystemGro):
        """
        Tests for some attributes and properties.
        """
        assert system_bmimbf4_gro.n_atoms == 1200

        test_mat = np.diag([4.26814, 4.26814, 4.26814])
        assert np.isclose(system_bmimbf4_gro.box_matrix, test_mat).all()

        test_line = 'ionic liquid'
        assert system_bmimbf4_gro.comment_line.strip() == test_line

        for index, mol_info_index in enumerate(system_bmimbf4_gro.molecules_info_ordered_all):
            if index < 300:
                assert mol_info_index == 0
            else:
                assert mol_info_index == 1

        test_mol_len = {
            ('BF4', 1): 0,
            ('BMIM', 3): 1
        }
        assert system_bmimbf4_gro.molecules_resname_len_index == test_mol_len

        test_comp = {
            'BF4': 300,
            'BMIM': 300
        }
        assert system_bmimbf4_gro.composition == test_comp
        assert 'Simulation system with:' in str(system_bmimbf4_gro)
        assert 'BMIM  : 300' in str(system_bmimbf4_gro)
        assert 'BF4   : 300' in str(system_bmimbf4_gro)

    def test_indexing(self, system_bmimbf4_gro: SystemGro):
        """
        Tests SystemGro accessing.
        """
        mol0 = system_bmimbf4_gro[0]
        mol1 = system_bmimbf4_gro[1]

        assert system_bmimbf4_gro[:2] == [mol0, mol1]
        assert system_bmimbf4_gro[-len(system_bmimbf4_gro)] == mol0
        assert system_bmimbf4_gro[-len(system_bmimbf4_gro):-len(system_bmimbf4_gro)+2] == [mol0, mol1]

        assert not system_bmimbf4_gro[200000:300000]
        assert not system_bmimbf4_gro[-300000:-200000]

        with pytest.raises(IndexError):
            _ = system_bmimbf4_gro[200000]
        with pytest.raises(IndexError):
            _ = system_bmimbf4_gro[-200000]

        with pytest.raises(TypeError):
            _ = system_bmimbf4_gro['test']  # type: ignore


class TestSystem:
    """
    Wraps all the tests for the System class.
    """
    def test_initialization(self):
        """
        Test System initialization with a system with molecules with multiple
        residues.
        """
        fgro = os.path.join(ACTUAL_PATH, '../../data/system_CG.gro')
        fgro = os.path.join(ACTUAL_PATH, '../../data/system_CG.gro')
        fitpSDS = os.path.join(ACTUAL_PATH, '../../data/SDS_AA.itp')
        fitpADN = os.path.join(ACTUAL_PATH, '../../data/DNA_CG.itp')

        with pytest.raises(IOError):
            _ = System(fgro, fitpSDS)

        empty = System(fgro)
        assert not empty.different_molecules
        with pytest.raises(IOError, match=r'The molecule.*'):
            empty.add_molecule_top(MoleculeTop(fitpSDS))

        bad_dna_top = MoleculeTop(fitpADN)

        for atom in bad_dna_top[98:105]:
            atom.resname = 'DA'
        for atom in bad_dna_top[105:112]:
            atom.resname = 'DG'

        with pytest.raises(IOError, match=r'The sequence.*'):
            empty.add_molecule_top(bad_dna_top)

        empty.add_molecule_top(MoleculeTop(fitpADN))

    def test_props(self, system: System):
        """
        Tests some system properties.
        """
        assert len(system) == 1068
        assert system[0].name == 'DNA'
        assert system[1].name == 'DPSM'
        assert system[-1].name == 'VTE'

        fgro = os.path.join(ACTUAL_PATH, '../../data/system_CG.gro')
        assert system.fgro == fgro

        with pytest.raises(IndexError):
            _ = system[-2000]
        with pytest.raises(IndexError):
            _ = system[2000]

        assert 'Simulation system with:' in str(system)
        assert 'DNA   : 1' in str(system)
        assert 'DPSM  : 67' in str(system)
        assert 'VTE   : 1000' in str(system)

        # Slice
        dna, dpsm = system[:2]
        assert dna.name == 'DNA'
        assert len(dna) == 232
        assert dpsm.name == 'DPSM'
        assert len(dpsm) == 9

        dna = system.different_molecules[0]
        assert dna == system[0]
        assert dna is not system[0]

        comps = system.composition
        assert comps['DNA'] == 1
        assert comps['DPSM'] == 67
        assert comps['VTE'] == 1000

    def test_indexing(self, system: System):
        """
        Tests System accessing.
        """
        mol0 = system[0]
        mol1 = system[1]

        assert system[:2] == [mol0, mol1]
        assert system[-len(system)] == mol0
        assert system[-len(system):-len(system)+2] == [mol0, mol1]

        assert not system[200000:300000]
        assert not system[-300000:-200000]

        with pytest.raises(IndexError):
            _ = system[200000]
        with pytest.raises(IndexError):
            _ = system[-200000]

        with pytest.raises(TypeError):
            _ = system['test']  # type: ignore
