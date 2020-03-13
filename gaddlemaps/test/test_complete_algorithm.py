'''
This module tests the whole Gaddle maps, starting with the creation of a system
to then calculate the alignment and finally map the system and write the final
.gro file.
'''

import os

import numpy as np
import pytest

from gaddlemaps import Manager
from gaddlemaps.components import System

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

sys_cg_gro = os.path.join(ACTUAL_PATH, '../data/sistema_CG.gro')

vte_cg_itp = os.path.join(ACTUAL_PATH, '../data/vitamin_E_CG.itp')
vte_cg_map_gro = os.path.join(ACTUAL_PATH, '../data/VTE_map.gro')

dna_cg_itp = os.path.join(ACTUAL_PATH, '../data/DNA_CG.itp')
dna_cg_map_gro = os.path.join(ACTUAL_PATH, '../data/DNA_map.gro')

vte_aa_itp = os.path.join(ACTUAL_PATH, '../data/VTE_AA.itp')
vte_aa_gro = os.path.join(ACTUAL_PATH, '../data/VTE_AA.gro')

dna_aa_itp = os.path.join(ACTUAL_PATH, '../data/DNA_AA.itp')
dna_aa_gro = os.path.join(ACTUAL_PATH, '../data/DNA_AA.gro')

bmim_aa_itp = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.itp')
bmim_aa_gro = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.gro')

bf4_aa_itp = os.path.join(ACTUAL_PATH, '../data/BF4_AA.itp')
bf4_aa_gro = os.path.join(ACTUAL_PATH, '../data/BF4_AA.gro')

bmim_cg_itp = os.path.join(ACTUAL_PATH, '../data/BMIM_CG.itp')
bf4_cg_itp = os.path.join(ACTUAL_PATH, '../data/BF4_CG.itp')

bmimbf4_sys = os.path.join(ACTUAL_PATH, '../data/system_bmimbf4_cg.gro')


def test_map_bmimbf4():
    sys = System(bmimbf4_sys, bmim_cg_itp, bf4_cg_itp)
    bmim = System(bmim_aa_gro, bmim_aa_itp)[0]
    bf4 = System(bf4_aa_gro, bf4_aa_itp)[0]
    man = Manager(sys)
    man.add_end_molecule(bmim)
    man.add_end_molecule(bf4)
    man.align_molecules()
    man.calculate_exchange_maps()

    fname = os.path.join(ACTUAL_PATH, '../data/sistema_CG_mapeado_test.gro')
    try:
        man.extrapolate_system(fname)
    except:
        import sys
        print('Error: ' + sys.exc_info()[0])
        try:
            os.remove(fname)
        except FileNotFoundError:
            pass
        raise

    fgro = fname
    sys_end = System(fgro, bmim_aa_itp, bf4_aa_itp)
    os.remove(fname)

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


@pytest.mark.longrun
def test_map_system():
    """
    Tests the whole maping process.
    """

    sys = System(sys_cg_gro, vte_cg_itp, dna_cg_itp)

    man = Manager(sys)

    svte = System(vte_aa_gro, vte_aa_itp)[0]
    sdna = System(dna_aa_gro, dna_aa_itp)[0]

    man.add_end_molecules(svte, sdna)

    start_dna = System(dna_cg_map_gro, dna_cg_itp)[0]
    start_vte = System(vte_cg_map_gro, vte_cg_itp)[0]

    # Tolerance in comparison
    poss = sys[0].atoms_positions
    tol = 0.1 * np.linalg.norm(np.max(poss, 0) + np.min(poss, 0))
    # Test aligning and not
    for i in range(2):
        if i:
            man.align_molecules()
        else:
            man.molecule_correspondence['VTE'].start = start_vte
            man.molecule_correspondence['DNA'].start = start_dna

        man.calculate_exchange_maps()

        fname = os.path.join(
            ACTUAL_PATH, '../data/sistema_CG_mapeado_test.gro')
        try:
            man.extrapolate_system(fname)
        except:
            import sys
            print('Error: ' + sys.exc_info()[0])
            try:
                os.remove(fname)
            except FileNotFoundError:
                pass
            raise

        fgro = fname
        sys_end = System(fgro, vte_aa_itp, dna_aa_itp)
        os.remove(fname)

        assert len(sys) == len(sys_end)
        assert sys.composition == sys_end.composition
        # The geometric center of each molecule should be more or less equal
        for mol1, mol2 in zip(sys, sys_end):
            assert np.allclose(mol1.geometric_center, mol2.geometric_center,
                               rtol=tol)
