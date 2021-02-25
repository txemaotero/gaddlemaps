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
Test for the _alignment submodule.
'''

import json
import os
from pathlib import Path

import numpy as np
import pytest

from gaddlemaps import (Alignment, ExchangeMap, guess_protein_restrains,
                        guess_residue_restrains, remove_hydrogens)
from gaddlemaps.components import (AtomGro, Molecule, MoleculeTop, Residue,
                                   System, SystemGro)
from gaddlemaps.parsers import GroFile

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def molecule_aa() -> Molecule:
    """
    Molecule instance of curcumine all atom.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/CUR_AA.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/CUR_AA.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def molecule_cg() -> Molecule:
    """
    Molecule instance of curcumine coarse-grained.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/CUR_map.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/CUR_CG.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def protein_aa() -> Molecule:
    """
    Molecule instance of a protein all-atom.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/Protein_AA.gro')
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/Protein_AA.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def protein_cg() -> Molecule:
    """
    Molecule instance of a protein coarse-grained.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/Protein_CG.gro')
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/Protein_CG.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def vte_aa() -> Molecule:
    """
    Molecule instance of vitamin E all-atom.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_AA.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_AA.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def vte_cg() -> Molecule:
    """
    Molecule instance of vitamin E coarse-grained.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/vitamin_E_CG.itp')
    return System(fgro, fitp)[0]


@pytest.fixture
def vte_map_cg() -> Molecule:
    """
    Molecule instance of vitamin E coarse-grained coming from the map
    process.
    """
    fgro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/VTE_map.gro')
    with GroFile(fgro) as _file:
        atoms = [AtomGro(line) for line in _file]
    mgro = Residue(atoms)
    fitp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/vitamin_E_CG.itp')
    mitp = MoleculeTop(fitp)
    return Molecule(mitp, [mgro])


@pytest.fixture
def alignment(molecule_aa: Molecule, molecule_cg: Molecule) -> Alignment:
    """
    Alignment instance with the curcumine molecule.
    """
    return Alignment(start=molecule_cg, end=molecule_aa)


def test_remove_hydrogens(vte_aa: Molecule):
    """
    Test for the remove_hydrogens function.
    """
    restr = [(0, 0), (19, 50)]
    pos, restr = remove_hydrogens(vte_aa, restr)
    restr_test = [(0, 0), (18, 50)]
    assert restr_test == restr


class TestAlignment:
    """
    Wraps all the tests for the Alignment class
    """
    def test_initialization(self, molecule_aa: Molecule, molecule_cg: Molecule):
        """
        Test the alignment class initialization.
        """
        with pytest.raises(TypeError):
            ali = Alignment(start=1)  # type: ignore
        with pytest.raises(TypeError):
            ali = Alignment(end=1)  # type: ignore
        ali = Alignment()
        assert ali.start is None
        assert ali.end is None
        # Some raises
        with pytest.raises(ValueError):
            ali.align_molecules()
        with pytest.raises(ValueError):
            ali.init_exchange_map()
        with pytest.raises(ValueError):
            ali.write_comparative_gro()
        ali.start = molecule_aa
        assert ali.exchange_map is None
        assert ali.start == molecule_aa
        ali.end = molecule_cg
        assert ali.end == molecule_cg
        # Reset
        ali.end = molecule_cg

        with pytest.raises(ValueError):
            ali.start = molecule_cg
        with pytest.raises(ValueError):
            ali.end = molecule_aa

    def test_align(self, alignment: Alignment):
        """
        Test that the alignment engine runs with are already aligned
        molecules.
        """
        # Simulate few steps
        alignment.STEPS_FACTOR = 10
        restr = [(0, 0), (0, 1)]
        def_type = (0, 1, 2)
        ign_h = False
        # Uncomment this to test a short simulation
        alignment.align_molecules(restr, def_type, ign_h)

    def test_exchange_map(self, alignment: Alignment):
        """
        Test the exchange map initialization from the alignment.
        """
        alignment.init_exchange_map(scale_factor=1)
        assert isinstance(alignment.exchange_map, ExchangeMap)

    def test_write_comparative_gro(self, alignment: Alignment, tmp_path: Path):
        """
        Tests the write_comparative_gro method.
        """
        subdir = tmp_path / "alignment_test"
        subdir.mkdir()
        tname = str(subdir / "compare.gro")

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

    def test_redefine_start(self, vte_aa: Molecule, vte_cg: Molecule,
                            vte_map_cg: Molecule):
        """
        Test the setter of the start molecule in alignment.
        """
        ali = Alignment(start=vte_cg)
        ali.end = vte_aa
        # Change start
        ali.start = vte_map_cg

        for at1, at2 in zip(ali.start, vte_map_cg):
            assert np.array_equal(at1.position, at2.position)
        for at1, at2 in zip(ali.end, vte_aa):
            assert np.array_equal(at1.position, at2.position)

    def test_exchange_map_vte(self, vte_aa: Molecule, vte_cg: Molecule,
                              vte_map_cg: Molecule):
        """
        Test the setter of the end molecule in alignment.
        """
        ali = Alignment(start=vte_map_cg, end=vte_aa)

        # For typing
        assert isinstance(ali.end, Molecule)
        assert isinstance(ali.start, Molecule)

        # Check that mols are well near.
        assert np.allclose(ali.start.geometric_center,
                           ali.end.geometric_center, rtol=0.5, atol=0.5)

        ali.init_exchange_map(scale_factor=1)
        assert isinstance(ali.exchange_map, ExchangeMap)

        new = ali.exchange_map(vte_cg)
        assert np.allclose(new.geometric_center,
                           vte_cg.geometric_center, rtol=0.5, atol=0.5)

    def test_different_molecules(self, protein_cg: Molecule,
                                 molecule_aa: Molecule):
        """
        Tests alignment with different molecules.
        """
        ali = Alignment(protein_cg, molecule_aa)
        with pytest.raises(IOError):
            ali.align_molecules(auto_guess_protein_restrictions=True)
        ali = Alignment(molecule_aa, protein_cg)
        ali.align_molecules()


def test_gess_residue_restrains(vte_cg: Molecule, vte_aa: Molecule):
    """
    Test gess_residue_restrains function with simple molecule.
    """
    restr = guess_residue_restrains(vte_cg, vte_aa)

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
