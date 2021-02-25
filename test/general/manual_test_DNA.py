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
"""
This script contains a test for the mapping of the DNA to run manually and
check the results visually.
"""


import os
import time

from gaddlemaps import Alignment
from gaddlemaps.components import System


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


dna_aa_itp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_AA.itp')
dna_aa_gro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_AA.gro')


sys_cg_gro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
sys_cg_gro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG.gro')
dna_cg_itp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/DNA_CG.itp')
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

