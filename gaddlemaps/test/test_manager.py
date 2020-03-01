#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
This module contains test for _manager module.

'''

import os
import json
import pytest
import numpy as np

from gaddlemaps import (Manager, Alignment)
from gaddlemaps.components import (System, Molecule,
                                   SystemGro, MoleculeItp)


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def system():
    fgro = os.path.join(ACTUAL_PATH, '../data/sistema_CG.gro')
    fitpDNA = os.path.join(ACTUAL_PATH, '../data/DNA_CG.itp')
    fitpVTE = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
    return System(fgro, fitpDNA, fitpVTE)


@pytest.fixture
def manager(system):
    return Manager(system)


@pytest.fixture
def manager_added(manager, molecule_DNA_AA, molecule_VTE_AA):
    manager.add_end_molecules(molecule_DNA_AA, molecule_VTE_AA)
    return manager


@pytest.fixture
def molecule_VTE_AA():
    fitpVTE = os.path.join(ACTUAL_PATH, '../data/VTE_AA.itp')
    fgroVTE = os.path.join(ACTUAL_PATH, '../data/VTE_AA.gro')
    VTEgro = SystemGro(fgroVTE)[0]
    VTEitp = MoleculeItp(fitpVTE)
    return Molecule(VTEgro, VTEitp)


@pytest.fixture
def molecule_DNA_AA():
    fitpDNA = os.path.join(ACTUAL_PATH, '../data/DNA_AA.itp')
    fgroDNA = os.path.join(ACTUAL_PATH, '../data/DNA_AA.gro')
    return System(fgroDNA, fitpDNA)[0]


@pytest.fixture
def molecule_VTE_map():
    fitpVTE = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
    fgroVTE = os.path.join(ACTUAL_PATH, '../data/VTE_map.gro')
    VTEgro = SystemGro(fgroVTE)[0]
    VTEitp = MoleculeItp(fitpVTE)
    return Molecule(VTEgro, VTEitp)


@pytest.fixture
def molecule_DNA_map():
    fitpDNA = os.path.join(ACTUAL_PATH, '../data/DNA_CG.itp')
    fgroDNA = os.path.join(ACTUAL_PATH, '../data/DNA_map.gro')
    return System(fgroDNA, fitpDNA)[0]


@pytest.fixture
def molecule_popc_AA():
    fitppopc = os.path.join(ACTUAL_PATH, '../data/popc-AA.itp')
    fgropopc = os.path.join(ACTUAL_PATH, '../data/popc-AA.gro')
    popcgro = SystemGro(fgropopc)[0]
    popcitp = MoleculeItp(fitppopc)
    return Molecule(popcgro, popcitp)


class TestManager(object):
    """
    Test for Manager class.

    """

    def test_init_from_files(self):
        fgro = os.path.join(ACTUAL_PATH, '../data/sistema_CG.gro')
        fitpDNA = os.path.join(ACTUAL_PATH, '../data/DNA_CG.itp')
        fitpVTE = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
        man = Manager.from_files(fgro, fitpVTE, fitpDNA)
        assert isinstance(man.system, System)
        assert len(man.molecule_correspondence) == 2
        assert 'DNA' in man.molecule_correspondence
        assert 'VTE' in man.molecule_correspondence
        assert isinstance(man.molecule_correspondence['DNA'].start, Molecule)
        assert isinstance(man.molecule_correspondence['VTE'].start, Molecule)

    def test_attrs(self, manager):
        man = manager
        assert isinstance(man.system, System)
        assert len(man.molecule_correspondence) == 2
        assert 'DNA' in man.molecule_correspondence
        assert 'VTE' in man.molecule_correspondence
        assert isinstance(man.molecule_correspondence['DNA'].start, Molecule)
        assert isinstance(man.molecule_correspondence['VTE'].start, Molecule)

    def test_add_mols(self, manager, molecule_DNA_AA, molecule_popc_AA,
                      molecule_VTE_AA):
        with pytest.raises(TypeError):
            manager.add_end_molecule(4)
        with pytest.raises(KeyError):
            manager.add_end_molecule(molecule_popc_AA)
        manager.add_end_molecules(molecule_DNA_AA, molecule_VTE_AA)
        assert manager.molecule_correspondence['DNA'].end == molecule_DNA_AA
        assert manager.molecule_correspondence['VTE'].end == molecule_VTE_AA

    def test_manager_added(self, manager_added):
        comp = manager_added.complete_correspondence
        assert 'DNA' in comp
        assert 'VTE' in comp

    def test_parse_restrictions(self, manager_added):
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
            rest = manager.parse_restrictions({'DNA': (1, 2)})
        with pytest.raises(ValueError):
            rest = manager.parse_restrictions({'DNA': [(1, 2, 3), ]})
        # Not changes
        rest_comp = {
            'DNA': [(0, 1), ],
            'VTE': [(0, 1), ],
        }
        rest = {
            'DNA': [(1, 2), ],
            'VTE': [(1, 2), ],
        }
        rest = manager.parse_restrictions(rest)
        assert rest == rest_comp
        # Complete
        rest_comp = {
            'DNA': [(0, 1), ],
            'VTE': None,
        }
        rest = {
            'DNA': [(1, 2), ],
        }
        rest = manager.parse_restrictions(rest)
        assert rest == rest_comp
        # Protein
        rest = manager.parse_restrictions(None, guess_proteins=True)
        assert rest['VTE'] == None
        assert len(rest['DNA']) == 1142

    def test_parse_deformations(self, manager_added):
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
            _ = manager._parse_deformations({'DNA': 1})
        with pytest.raises(ValueError):
            _ = manager._parse_deformations({'DNA': (1, 2, 3, 4, 4)})
        # Not changes
        def_comp = {
            'DNA': (1, 2),
            'VTE': (1, 2),
        }
        defo = manager._parse_deformations(def_comp)
        assert defo == def_comp
        # Complete
        def_comp = {
            'DNA': (1, 2),
            'VTE': None,
        }
        defo = {
            'DNA': (1, 2),
        }
        defo = manager._parse_deformations(defo)
        assert defo == def_comp

    def test_parse_ignore_hydrogens(self, manager_added):
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
            _ = manager._parse_ignore_hydrogens({'DNA': 1})
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

    def test_info(self, manager_added, system, molecule_DNA_AA, molecule_VTE_AA):
        info = manager_added.info
        mols_per_name = system.different_molecules
        end_name = {'DNA': molecule_DNA_AA, 'VTE': molecule_VTE_AA}
        mol_corr_info = {mol.name: Alignment(mol, end_name[mol.name]).info
                         for mol in mols_per_name}
        info_test = {
            'system': system.info,
            'molecule_correspondence': mol_corr_info,
        }
        assert info_test == info
        new_man = Manager.from_info(info)
        assert len(new_man.system) == len(manager_added.system)
        assert new_man.molecule_correspondence.keys(
        ) == manager_added.molecule_correspondence.keys()
        for (n1, ali1), (n2, ali2) in zip(new_man.molecule_correspondence.items(),
                                          manager_added.molecule_correspondence.items()):
            assert n1 == n2
            assert ali1.start == ali2.start
            assert ali1.end == ali2.end
        # Load from json
        # json test
        text = json.dumps(info, indent=2)
        parse_text = json.loads(text)
        new_man = Manager.from_info(parse_text)
        assert len(new_man.system) == len(manager_added.system)
        assert sorted(new_man.molecule_correspondence.keys()) == sorted(
            manager_added.molecule_correspondence.keys())
        for (n1, ali1), (n2, ali2) in zip(sorted(new_man.molecule_correspondence.items()),
                                          sorted(manager_added.molecule_correspondence.items())):
            assert n1 == n2
            assert ali1.start == ali2.start
            assert ali1.end == ali2.end

    def test_extrapolate(self, manager, molecule_DNA_map, molecule_VTE_map,
                         system, molecule_VTE_AA, molecule_DNA_AA):
        fname = os.path.join(
            ACTUAL_PATH, '../data/sistema_CG_mapeado_test.gro')
        with pytest.raises(SystemError):
            _ = manager.extrapolate_system(fname)

        manager.molecule_correspondence['VTE'].start = molecule_VTE_map
        manager.molecule_correspondence['DNA'].start = molecule_DNA_map

        manager.add_end_molecules(molecule_DNA_AA, molecule_VTE_AA)

        tot = manager.complete_correspondence
        assert 'DNA' in tot
        assert 'VTE' in tot

        manager.calculate_exchange_maps()

        try:
            manager.extrapolate_system(fname)
        except:
            import sys
            print('Error: ' + sys.exc_info()[0])
            try:
                os.remove(fname)
            except FileNotFoundError:
                pass
            raise

        sys = system

        fgro = fname
        fitpDNA = os.path.join(ACTUAL_PATH, '../data/DNA_AA.itp')
        fitpVTE = os.path.join(ACTUAL_PATH, '../data/VTE_AA.itp')
        sys_end = System(fgro, fitpDNA, fitpVTE)
        os.remove(fname)

        poss = sys[0].atoms_positions
        tol = 0.1 * np.linalg.norm(np.max(poss, 0) + np.min(poss, 0))

        assert len(sys) == len(sys_end)
        assert sys.composition == sys_end.composition
        # The geometric center of each molecule should be more or less equal
        for mol1, mol2 in zip(sys, sys_end):
            assert np.allclose(mol1.geometric_center, mol2.geometric_center,
                               rtol=tol, atol=tol)
