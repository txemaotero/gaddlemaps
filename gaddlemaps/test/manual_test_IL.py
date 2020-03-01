"""
This script contains a test for the mapping of the DNA to run manually and
check the results visually.
"""


import os
import time

from gaddlemaps import Manager
from gaddlemaps import Alignment
from gaddlemaps.components import System


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

bmim_aa_itp = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.itp')
bmim_aa_gro = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.gro')

bf4_aa_itp = os.path.join(ACTUAL_PATH, '../data/BF4_AA.itp')
bf4_aa_gro = os.path.join(ACTUAL_PATH, '../data/BF4_AA.gro')

bmim_cg_itp = os.path.join(ACTUAL_PATH, '../data/BMIM_CG.itp')
bf4_cg_itp = os.path.join(ACTUAL_PATH, '../data/BF4_CG.itp')

bmimbf4_sys = os.path.join(ACTUAL_PATH, '../data/system_bmimbf4_cg.gro')


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
