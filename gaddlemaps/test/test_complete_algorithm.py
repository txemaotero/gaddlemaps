'''
This module tests the whole Gaddle maps, starting with the creation of a system
to then calculate the alignment and finally map the system and write the final
.gro file.
'''

import os
from pathlib import Path

import numpy as np
import pytest

from gaddlemaps import Manager
from gaddlemaps.components import System

# Define some file names
ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

SYS_CG_GRO = os.path.join(ACTUAL_PATH, '../data/sistema_CG.gro')

VTE_CG_ITP = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
VTE_CG_MAP_GRO = os.path.join(ACTUAL_PATH, '../data/VTE_map.gro')

DNA_CG_ITP = os.path.join(ACTUAL_PATH, '../data/DNA_CG.itp')
DNA_CG_MAP_GRO = os.path.join(ACTUAL_PATH, '../data/DNA_map.gro')

VTE_AA_ITP = os.path.join(ACTUAL_PATH, '../data/VTE_AA.itp')
VTE_AA_GRO = os.path.join(ACTUAL_PATH, '../data/VTE_AA.gro')

DNA_AA_ITP = os.path.join(ACTUAL_PATH, '../data/DNA_AA.itp')
DNA_AA_GRO = os.path.join(ACTUAL_PATH, '../data/DNA_AA.gro')

BMIM_AA_ITP = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.itp')
BMIM_AA_GRO = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.gro')

BF4_AA_ITP = os.path.join(ACTUAL_PATH, '../data/BF4_AA.itp')
BF4_AA_GRO = os.path.join(ACTUAL_PATH, '../data/BF4_AA.gro')

BMIM_CG_ITP = os.path.join(ACTUAL_PATH, '../data/BMIM_CG.itp')
BF4_CG_ITP = os.path.join(ACTUAL_PATH, '../data/BF4_CG.itp')

BMIMBF4_SYS = os.path.join(ACTUAL_PATH, '../data/system_bmimbf4_cg.gro')


def test_map_bmimbf4(tmp_path: Path):
    """
    Tests the map process for the BMIM BF4 system. 
    """
    sys = System(BMIMBF4_SYS, BMIM_CG_ITP, BF4_CG_ITP)
    bmim = System(BMIM_AA_GRO, BMIM_AA_ITP)[0]
    bf4 = System(BF4_AA_GRO, BF4_AA_ITP)[0]
    man = Manager(sys)
    man.add_end_molecule(bmim)
    man.add_end_molecule(bf4)
    man.align_molecules()
    man.calculate_exchange_maps()

    subdir = tmp_path / "complete_alg_test"
    subdir.mkdir()
    fname = str(subdir / "system_CG_mapped_test.gro")

    man.extrapolate_system(fname)

    sys_end = System(fname, BMIM_AA_ITP, BF4_AA_ITP)

    # Tolerance in comparison
    poss = sys[0].atoms_positions
    tol = 0.1 * np.linalg.norm(np.max(poss, 0) + np.min(poss, 0))

    assert len(sys) == len(sys_end)
    assert sys.composition == sys_end.composition
    # The geometric center of each molecule should be more or less equal
    for mol1, mol2 in zip(sys, sys_end):
        geom1 = mol1.geometric_center
        geom2 = mol2.geometric_center
        assert np.allclose(geom1, geom2, rtol=tol, atol=tol)


def test_map_system(tmp_path: Path):
    """
    Tests the whole mapping process.
    """
    sys = System(SYS_CG_GRO, VTE_CG_ITP, DNA_CG_ITP)

    man = Manager(sys)

    svte = System(VTE_AA_GRO, VTE_AA_ITP)[0]
    sdna = System(DNA_AA_GRO, DNA_AA_ITP)[0]

    man.add_end_molecules(svte, sdna)

    start_dna = System(DNA_CG_MAP_GRO, DNA_CG_ITP)[0]
    start_vte = System(VTE_CG_MAP_GRO, VTE_CG_ITP)[0]

    # Tolerance in comparison
    poss = sys[0].atoms_positions
    tol = 0.1 * np.linalg.norm(np.max(poss, 0) + np.min(poss, 0))

    subdir = tmp_path / "complete_alg_test"
    subdir.mkdir()
    fname = str(subdir / "system_CG_mapped_test.gro")

    # Test aligning and not
    for i in range(2):
        if i:
            man.align_molecules()
        else:
            man.molecule_correspondence['VTE'].start = start_vte
            man.molecule_correspondence['DNA'].start = start_dna

        man.calculate_exchange_maps()

        man.extrapolate_system(fname)

        sys_end = System(fname, VTE_AA_ITP, DNA_AA_ITP)

        assert len(sys) == len(sys_end)
        assert sys.composition == sys_end.composition
        # The geometric center of each molecule should be more or less equal
        for mol1, mol2 in zip(sys, sys_end):
            assert np.allclose(mol1.geometric_center, mol2.geometric_center,
                               rtol=tol)
