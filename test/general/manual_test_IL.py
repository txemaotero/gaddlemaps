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

from gaddlemaps import Manager
from gaddlemaps import Alignment
from gaddlemaps.components import System


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

bmim_aa_itp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/BMIM_AA.itp')
bmim_aa_gro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/BMIM_AA.gro')

bf4_aa_itp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/BF4_AA.itp')
bf4_aa_gro = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/BF4_AA.gro')

bmim_cg_itp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/BMIM_CG.itp')
bf4_cg_itp = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/BF4_CG.itp')

bmimbf4_sys = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_bmimbf4_cg.gro')


system = System(bmimbf4_sys, bmim_cg_itp, bf4_cg_itp)
bmim = System(bmim_aa_gro, bmim_aa_itp)[0]
bf4 = System(bf4_aa_gro, bf4_aa_itp)[0]
man = Manager(system)
man.add_end_molecule(bmim)
man.add_end_molecule(bf4)
man.align_molecules()
man.calculate_exchange_maps()

fname = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG_mapeado_test.gro')
fname = os.path.join(ACTUAL_PATH, '../../gaddlemaps/data/system_CG_mapeado_test.gro')
try:
    man.extrapolate_system(fname)
except:
    import sys
    print(f'Error: {sys.exc_info()[0]}')
    try:
        os.remove(fname)
    except FileNotFoundError:
        pass
    raise
