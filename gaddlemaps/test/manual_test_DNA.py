"""
This script contains a test for the mapping of the DNA to run manually and
check the results visually.
"""


import os
import time

from gaddlemaps import Alignment
from gaddlemaps.components import System


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


dna_aa_itp = os.path.join(ACTUAL_PATH, '../data/DNA_AA.itp')
dna_aa_gro = os.path.join(ACTUAL_PATH, '../data/DNA_AA.gro')


sys_cg_gro = os.path.join(ACTUAL_PATH, '../data/sistema_CG.gro')
dna_cg_itp = os.path.join(ACTUAL_PATH, '../data/DNA_CG.itp')
sys = System(sys_cg_gro, dna_cg_itp)

Protein_cg = sys[0]
Protein_aa = System(dna_aa_gro, dna_aa_itp)[0]

ali = Alignment(Protein_cg, Protein_aa)

start = time.time()
ali.align_molecules()
end = time.time()

time_i = end - start
print(f'Guessed restrains: {time_i}')

ali.write_comparative_gro('check_dna.gro')

