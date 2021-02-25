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
This module contains test for _manager module.
'''

import json
import os
from pathlib import Path

import numpy as np
import pytest

from gaddlemaps import Alignment, Manager
from gaddlemaps.components import Molecule, MoleculeTop, System, SystemGro

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def system() -> System:
    """
    System instance with DNA and E vitamin E molecules.
    """

    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
    fitpDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_CG.itp')
    fitpVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/vitamin_E_CG.itp')
    return System(fgro, fitpDNA, fitpVTE)


@pytest.fixture
def molecule_VTE_AA() -> Molecule:
    """
    Molecule instance of E vitamin in all-atom resolution.
    """
    fitpVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_AA.itp')
    fgroVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_AA.gro')
    VTEgro = SystemGro(fgroVTE)[0]
    VTEitp = MoleculeTop(fitpVTE)
    return Molecule(VTEitp, [VTEgro])


@pytest.fixture
def molecule_VTE_map() -> Molecule:
    """
    Molecule instance of E vitamin in coarse-grained resolution after being
    mapped.
    """
    fitpVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/vitamin_E_CG.itp')
    fgroVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_map.gro')
    VTEgro = SystemGro(fgroVTE)[0]
    VTEitp = MoleculeTop(fitpVTE)
    return Molecule(VTEitp, [VTEgro])


@pytest.fixture
def molecule_DNA_AA() -> Molecule:
    """
    Molecule instance of DNA in all-atom resolution.
    """
    fitpDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_AA.itp')
    fgroDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_AA.gro')
    return System(fgroDNA, fitpDNA)[0]


@pytest.fixture
def molecule_DNA_map() -> Molecule:
    """
    Molecule instance of DNA in coarse-grained resolution after being mapped.
    """
    fitpDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_CG.itp')
    fgroDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_map.gro')
    return System(fgroDNA, fitpDNA)[0]


@pytest.fixture
def molecule_popc_AA() -> Molecule:
    """
    Molecule instance of POPC in all-atom resolution.
    """
    fitppopc = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/popc-AA.itp')
    fgropopc = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/popc-AA.gro')
    popcgro = SystemGro(fgropopc)[0]
    popcitp = MoleculeTop(fitppopc)
    return Molecule(popcitp, [popcgro])


@pytest.fixture
def manager(system: System) -> Manager:
    """
    Manager instance initialized with the system with DNA and E vitamin but
    without molecule assignment in the final resolution.
    """
    return Manager(system)


@pytest.fixture
def manager_added(manager: Manager, molecule_DNA_AA: Molecule,
                  molecule_VTE_AA: Molecule) -> Manager:
    """
    Manager instance initialized with the system with DNA and E vitamin
    with molecules in the final resolution assigned.
    """
    manager.add_end_molecules(molecule_DNA_AA, molecule_VTE_AA)
    return manager


class TestManager:
    """
    Wraps all tests for Manager class.
    """
    def test_init_from_files(self):
        """
        Tests the Manager initialization from corresponding files.
        """
        fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
        fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
        fitpDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_CG.itp')
        fitpVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/vitamin_E_CG.itp')
        man = Manager.from_files(fgro, fitpVTE, fitpDNA)
        assert isinstance(man.system, System)
        assert len(man.molecule_correspondence) == 2
        assert 'DNA' in man.molecule_correspondence
        assert 'VTE' in man.molecule_correspondence
        assert isinstance(man.molecule_correspondence['DNA'].start, Molecule)
        assert isinstance(man.molecule_correspondence['VTE'].start, Molecule)

    def test_attributes(self, manager: Manager):
        """
        Tests some attributes with an initialized Manager.
        """
        man = manager
        assert isinstance(man.system, System)
        assert len(man.molecule_correspondence) == 2
        assert 'DNA' in man.molecule_correspondence
        assert 'VTE' in man.molecule_correspondence
        assert isinstance(man.molecule_correspondence['DNA'].start, Molecule)
        assert isinstance(man.molecule_correspondence['VTE'].start, Molecule)

    def test_add_mols(self, manager: Manager, molecule_DNA_AA: Molecule,
                      molecule_popc_AA: Molecule, molecule_VTE_AA: Molecule):
        """
        Tests the final resolution Molecule assignment.
        """
        with pytest.raises(TypeError):
            manager.add_end_molecule(4)  # type: ignore
        with pytest.raises(KeyError):
            manager.add_end_molecule(molecule_popc_AA)
        manager.add_end_molecules(molecule_DNA_AA, molecule_VTE_AA)
        assert manager.molecule_correspondence['DNA'].end == molecule_DNA_AA
        assert manager.molecule_correspondence['VTE'].end == molecule_VTE_AA

    def test_manager_added(self, manager_added: Manager):
        """
        Tests manager complete_correspondence attribute.
        """
        comp = manager_added.complete_correspondence
        assert 'DNA' in comp
        assert 'VTE' in comp

    def test_parse_restrictions(self, manager_added: Manager):
        """
        Tests pares_restriction methods.
        """
        manager = manager_added

        # If its None
        rest_comp = {
            'DNA': None,
            'VTE': None,
        }
        rest = manager.parse_restrictions(None)
        assert rest == rest_comp

        # Wrong name
        with pytest.raises(KeyError):
            rest = manager.parse_restrictions({'test': [(1, 2), ]})

        # Wrong format
        with pytest.raises(ValueError):
            rest = manager.parse_restrictions({'DNA': (1, 2)})  # type: ignore
        with pytest.raises(ValueError):
            rest = manager.parse_restrictions({'DNA': [(1, 2, 3), ]})  # type: ignore

        # No changes
        rest_comp1 = {
            'DNA': [(0, 1), ],
            'VTE': [(0, 1), ],
        }
        new_rest = manager.parse_restrictions(rest_comp1)  # type: ignore
        assert new_rest == rest_comp1

        # Complete
        rest_comp2 = {
            'DNA': [(0, 1), ],
            'VTE': None,
        }
        rest = {
            'DNA': [(0, 1), ],
        }
        rest = manager.parse_restrictions(rest)
        assert rest == rest_comp2
        # Protein
        rest = manager.parse_restrictions(None, guess_proteins=True)
        assert rest['VTE'] == None
        assert len(rest['DNA']) == 1142

    def test_parse_deformations(self, manager_added: Manager):
        """
        Tests the _parse_deformations method
        """
        manager = manager_added
        # If its None
        def_comp = {
            'DNA': None,
            'VTE': None,
        }
        defo = manager._parse_deformations(None)
        assert defo == def_comp
        # Wrong name
        with pytest.raises(KeyError):
            _ = manager._parse_deformations({'test': (2, 0)})
        # Wrong format
        with pytest.raises(ValueError):
            _ = manager._parse_deformations({'DNA': 1})  # type: ignore
        with pytest.raises(ValueError):
            _ = manager._parse_deformations({'DNA': (1, 2, 3, 4, 4)})
        # Not changes
        def_comp1 = {
            'DNA': (1, 2),
            'VTE': (1, 2),
        }
        defo = manager._parse_deformations(def_comp1)  # type: ignore
        assert defo == def_comp1
        # Complete
        def_comp2 = {
            'DNA': (1, 2),
            'VTE': None,
        }
        defo = {
            'DNA': (1, 2),
        }
        defo = manager._parse_deformations(defo)
        assert defo == def_comp2

    def test_parse_ignore_hydrogens(self, manager_added: Manager):
        """
        Tests _parse_ignore_hydrogens method.
        """
        manager = manager_added
        # If its None
        ign_comp = {
            'DNA': True,
            'VTE': True,
        }
        ign = manager._parse_ignore_hydrogens(None)
        assert ign == ign_comp
        with pytest.raises(KeyError):
            _ = manager._parse_ignore_hydrogens({'test': True})
        # Wrong format
        with pytest.raises(ValueError):
            _ = manager._parse_ignore_hydrogens({'DNA': 1})  # type: ignore
        # Complete
        ign_comp = {
            'DNA': False,
            'VTE': True,
        }
        ign = {
            'DNA': False,
        }
        ign = manager._parse_ignore_hydrogens(ign)
        assert ign == ign_comp

    def test_extrapolate(self, manager: Manager, molecule_DNA_map: Molecule,
                         molecule_VTE_map: Molecule, system: System,
                         molecule_VTE_AA: Molecule, molecule_DNA_AA: Molecule,
                         tmp_path: Path):
        """
        Tests the extrapolation functionality.
        """
        subdir = tmp_path / "manager_test"
        subdir.mkdir()
        fname = str(subdir / "system_CG_mapped_test.gro")

        with pytest.raises(SystemError):
            _ = manager.extrapolate_system(fname)

        manager.molecule_correspondence['VTE'].start = molecule_VTE_map
        manager.molecule_correspondence['DNA'].start = molecule_DNA_map

        manager.add_end_molecules(molecule_DNA_AA, molecule_VTE_AA)

        tot = manager.complete_correspondence
        assert 'DNA' in tot
        assert 'VTE' in tot

        manager.calculate_exchange_maps()

        manager.extrapolate_system(fname)

        fgro = fname
        fitpDNA = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_AA.itp')
        fitpVTE = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_AA.itp')
        sys_end = System(fgro, fitpDNA, fitpVTE)

        poss = system[0].atoms_positions
        tol = 0.1 * np.linalg.norm(np.max(poss, 0) + np.min(poss, 0))

        assert len(system) == len(sys_end)
        assert system.composition == sys_end.composition
        # The geometric center of each molecule should be more or less equal
        for mol1, mol2 in zip(system, sys_end):
            assert np.allclose(mol1.geometric_center, mol2.geometric_center,
                               rtol=tol, atol=tol)
