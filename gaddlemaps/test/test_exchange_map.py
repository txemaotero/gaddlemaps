#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Test for exchange map module.

'''


import pytest
import os
import numpy as np
from gaddlemaps.components import (Molecule, MoleculeGro, MoleculeItp,
                                   AtomGro, System)
from gaddlemaps.parsers import GroFile
from gaddlemaps import ExchangeMap


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def molecule_aa():
    fgro = os.path.join(ACTUAL_PATH, '../data/CUR_AA.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = MoleculeGro(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/CUR_AA.itp')
    mitp = MoleculeItp(fitp)
    return Molecule(mgro, mitp)


@pytest.fixture
def molecule_cg():
    fgro = os.path.join(ACTUAL_PATH, '../data/CUR_map.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = MoleculeGro(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/CUR_CG.itp')
    mitp = MoleculeItp(fitp)
    return Molecule(mgro, mitp)


@pytest.fixture
def bf4_aa():
    fgro = os.path.join(ACTUAL_PATH, '../data/BF4_AA.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/BF4_AA.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def bf4_cg():
    fgro = os.path.join(ACTUAL_PATH, '../data/BF4_CG.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/BF4_CG.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def vte_aa():
    fgro = os.path.join(ACTUAL_PATH, '../data/VTE_AA.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = MoleculeGro(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/VTE_AA.itp')
    mitp = MoleculeItp(fitp)
    return Molecule(mgro, mitp)


@pytest.fixture
def vte_cg():
    fgro = os.path.join(ACTUAL_PATH, '../data/sistema_CG.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def vte_map_cg():
    fgro = os.path.join(ACTUAL_PATH, '../data/VTE_map.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = MoleculeGro(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
    mitp = MoleculeItp(fitp)
    return Molecule(mgro, mitp)


def test_exchange_map(molecule_cg, molecule_aa):
    """
    Test the ExchangeMap class.

    """
    emap = ExchangeMap(molecule_cg, molecule_aa, scale_factor=1)
    with pytest.raises(TypeError, match='.*Argument.*'):
        emap(2)
    with pytest.raises(TypeError, match='.*refmolecule.*'):
        emap(molecule_aa)
    map_mol = emap(molecule_cg)
    assert isinstance(map_mol, Molecule)
    assert len(map_mol) == 41
    assert map_mol == molecule_aa
    assert map_mol is not molecule_aa
    for at1, at2 in zip(map_mol, molecule_aa):
        assert np.allclose(at1.position, at2.position)
    # Reverse map
    emap = ExchangeMap(molecule_aa, molecule_cg, scale_factor=1)
    map_mol = emap(molecule_aa)
    assert isinstance(map_mol, Molecule)
    assert len(map_mol) == 8
    assert map_mol == molecule_cg
    assert map_mol is not molecule_cg
    for at1, at2 in zip(map_mol, molecule_cg):
        assert np.allclose(at1.position, at2.position)


def test_bf4_map(bf4_cg, bf4_aa):
    emap = ExchangeMap(bf4_cg, bf4_aa, scale_factor=1)
    with pytest.raises(TypeError, match='.*Argument.*'):
        emap(2)
    with pytest.raises(TypeError, match='.*refmolecule.*'):
        emap(bf4_aa)
    map_mol = emap(bf4_cg)
    assert isinstance(map_mol, Molecule)
    assert len(map_mol) == 5
    assert map_mol == bf4_aa
    assert map_mol is not bf4_aa
    assert np.allclose(map_mol.geometric_center, bf4_aa.geometric_center,
                       rtol=1.e-3, atol=1.e-3)
    # Reverse map
    emap = ExchangeMap(bf4_aa, bf4_cg, scale_factor=1)
    map_mol = emap(bf4_aa)
    assert isinstance(map_mol, Molecule)
    assert len(map_mol) == 1
    assert map_mol == bf4_cg
    assert map_mol is not bf4_cg
    assert np.allclose(map_mol.geometric_center, bf4_cg.geometric_center,
                       rtol=1.e-3, atol=1.e-3)


def test_vte_map(vte_aa, vte_cg, vte_map_cg):
    emap = ExchangeMap(vte_map_cg, vte_aa, scale_factor=1)
    # Check the copy of refmolecule
    start_positions = vte_aa.atoms_positions
    new = emap(vte_cg)
    end_positions = vte_aa.atoms_positions
    assert np.array_equal(start_positions, end_positions)
    assert new is not vte_aa
    assert len(new) == len(vte_aa)
    for at1, at2 in zip(new, vte_aa):
        assert not np.array_equal(at1.position, at2.position)
    for pos1, pos2 in zip(new.atoms_positions, vte_aa.atoms_positions):
        assert not np.array_equal(pos1, pos2)
    assert not np.array_equal(new.geometric_center, vte_aa.geometric_center)
    assert np.allclose(new.geometric_center, vte_cg.geometric_center,
                       rtol=0.5, atol=0.5)
