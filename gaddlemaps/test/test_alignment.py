'''
Test for the _alignment submodule.
'''

import json
import os

import numpy as np
import pytest

from gaddlemaps import (Alignment, ExchangeMap, guess_molecule_restrains,
                        guess_protein_restrains, remove_hydrogens)
from gaddlemaps.components import (AtomGro, Molecule, MoleculeTop, Residue,
                                   System, SystemGro)
from gaddlemaps.parsers import GroFile

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def molecule_aa():
    fgro = os.path.join(ACTUAL_PATH, '../data/CUR_AA.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/CUR_AA.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def molecule_cg():
    fgro = os.path.join(ACTUAL_PATH, '../data/CUR_map.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/CUR_CG.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def protein_aa():
    fgro = os.path.join(ACTUAL_PATH, '../data/Protein_AA.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/Protein_AA.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def protein_cg():
    fgro = os.path.join(ACTUAL_PATH, '../data/Protein_CG.gro')
    fitp = os.path.join(ACTUAL_PATH, '../data/Protein_CG.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def vte_aa():
    fgro = os.path.join(ACTUAL_PATH, '../data/VTE_AA.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/VTE_AA.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


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
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def alignment(molecule_aa, molecule_cg):
    return Alignment(start=molecule_cg, end=molecule_aa)


def test_remove_hydrogens(vte_aa):
    restr = [(0, 0), (19, 50)]
    pos, restr = remove_hydrogens(vte_aa, restr)
    restr_test = [(0, 0), (18, 50)]
    assert restr_test == restr


class TestAlignment(object):
    """
    Test for the Alignment class

    """

    def test_init(self, molecule_aa, molecule_cg):
        with pytest.raises(TypeError):
            ali = Alignment(start=1)
        with pytest.raises(TypeError):
            ali = Alignment(end=1)
        ali = Alignment()
        assert ali.start is None
        assert ali.end is None
        ali.start = molecule_aa
        assert ali.exchange_map is None
        assert ali.start == molecule_aa
        ali.end = molecule_cg
        with pytest.raises(ValueError):
            ali.start = molecule_cg
        with pytest.raises(ValueError):
            ali.end = molecule_aa

    def test_align(self, alignment):
        """
        The molecules are already aligned.
        """
        # Simulate few steps
        alignment.STEPS_FACTOR = 10
        restr = [(0, 0), (0, 1)]
        def_type = (0, 1, 2)
        ign_h = False
        # Uncomment this to test a short simulation
        # alignment.align_molecules(restr, def_type, ign_h)

    def test_exchange_map(self, alignment):
        alignment.init_exchange_map(scale_factor=1)
        assert isinstance(alignment.exchange_map, ExchangeMap)

    def test_compare_file(self, alignment):
        tname = 'compare.gro'
        alignment.write_comparative_gro(tname)
        # Load the sys
        s = SystemGro(tname)
        assert len(s) == 2
        start, end = s[:2]
        assert start.resname == 'START'
        assert end.resname == 'END'
        assert len(start) == 8
        assert len(end) == 41
        os.remove(tname)

    def test_redefine_start(self, vte_aa, vte_cg, vte_map_cg):
        ali = Alignment(start=vte_cg)
        ali.end = vte_aa
        # Change start
        ali.start = vte_map_cg

        for at1, at2 in zip(ali.start, vte_map_cg):
            assert np.array_equal(at1.position, at2.position)
        for at1, at2 in zip(ali.end, vte_aa):
            assert np.array_equal(at1.position, at2.position)

    def test_exchange_map_vte(self, vte_aa, vte_cg, vte_map_cg):
        ali = Alignment(start=vte_map_cg, end=vte_aa)

        # Check that mols are well near.
        assert np.allclose(ali.start.geometric_center,
                           ali.end.geometric_center, rtol=0.5, atol=0.5)

        ali.init_exchange_map(scale_factor=1)
        new = ali.exchange_map(vte_cg)
        assert np.allclose(new.geometric_center,
                           vte_cg.geometric_center, rtol=0.5, atol=0.5)

    def test_info(self, alignment, molecule_aa, molecule_cg):
        info = alignment.info
        info_test = {
            'start': molecule_cg.info.to_dict(),
            'end': molecule_aa.info.to_dict(),
            'exchange_map': None,
        }
        assert info == info_test
        alignment.init_exchange_map(0.6)
        info = alignment.info
        info_test['exchange_map'] = 0.6
        assert info == info_test
        new_alig = Alignment.from_info(info)
        assert alignment.start == new_alig.start
        assert alignment.end == new_alig.end
        assert alignment.exchange_map.scale_factor == alignment.exchange_map.scale_factor
        # json test
        text = json.dumps(info, indent=2)
        parse_text = json.loads(text)
        new_alig = Alignment.from_info(parse_text)
        assert alignment.start == new_alig.start
        assert alignment.end == new_alig.end
        assert alignment.exchange_map.scale_factor == alignment.exchange_map.scale_factor


def test_gess_molecule_restrains(vte_cg, vte_aa):
    """
    Test gess_molecule_restrains function with simple molecule.
    """
    restr = guess_molecule_restrains(vte_cg, vte_aa)

    restr_test = [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6),
                  (2, 7), (2, 8), (3, 9), (3, 10), (3, 11), (4, 12), (4, 13),
                  (4, 14), (4, 15), (5, 16), (5, 17), (5, 18), (6, 19),
                  (6, 20), (6, 21), (7, 22), (7, 23), (7, 24), (8, 25),
                  (8, 26), (8, 27), (9, 28), (9, 29), (9, 30), (9, 31)]
    assert restr_test == restr


def test_gess_protein_restrains(protein_cg, protein_aa):
    """
    Test gess_protein_restrains function with simple molecule.
    """
    restr = guess_protein_restrains(protein_cg, protein_aa)
    assert (0, 0) in restr
    assert (8, 47) in restr
    assert (23, 126) in restr
