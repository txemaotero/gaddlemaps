# -*- coding: utf-8 -*-
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
This module run molecular alignment to calculate the calculation time.
"""


import os
import time
import numpy as np

from gaddlemaps import Alignment
from gaddlemaps.components import System


ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]

Protein_cg_gro = os.path.join(ACTUAL_PATH, '../data/Protein_CG.gro')
Protein_cg_itp = os.path.join(ACTUAL_PATH, '../data/Protein_CG.itp')

Protein_aa_gro = os.path.join(ACTUAL_PATH, '../data/Protein_AA.gro')
Protein_aa_itp = os.path.join(ACTUAL_PATH, '../data/Protein_AA.itp')

n_test = 1
times = np.empty(n_test)

for i in range(n_test):
    Protein_cg = System(Protein_cg_gro, Protein_cg_itp)[0]
    Protein_aa = System(Protein_aa_gro, Protein_aa_itp)[0]

    ali = Alignment(Protein_cg, Protein_aa)
    start = time.time()
    ali.align_molecules()
    end = time.time()
    time_i = end - start
    times[i] = time_i
    print(f'Guessed restrains: {time_i}')

avg = np.mean(times)
std = np.std(times)
print(f'Average time {avg} +- {std}')


times = np.empty(n_test)

for i in range(n_test):
    Protein_cg = System(Protein_cg_gro, Protein_cg_itp)[0]
    Protein_aa = System(Protein_aa_gro, Protein_aa_itp)[0]

    ali = Alignment(Protein_cg, Protein_aa)
    start = time.time()
    ali.align_molecules(restrictions=[(0, 0)])
    end = time.time()
    time_i = end - start
    times[i] = time_i
    print(f'One restrain: {time_i}')

for i in range(n_test):
    Protein_cg = System(Protein_cg_gro, Protein_cg_itp)[0]
    Protein_aa = System(Protein_aa_gro, Protein_aa_itp)[0]

    ali = Alignment(Protein_cg, Protein_aa)
    start = time.time()
    ali.align_molecules(restrictions=[])
    end = time.time()
    time_i = end - start
    times[i] = time_i
    print(f'No restrain: {time_i}')

avg = np.mean(times)
std = np.std(times)
print(f'Average time {avg} +- {std}')
