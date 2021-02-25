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
'''
gaddlemaps --- Change molecules in simulations with python
==========================================================

This python package provides an implementation of the GADDLE-Maps (General
Algorithm for Discrete Object Deformations Based on Local Exchange Maps) that
allows to change molecules in a molecular dynamics simulations. For example,
this tool is very handy to back-map a coarse grained simulation to atomistic
forcefields and also in the other way around.

TODO: Complete this description
'''



import os
from typing import Optional, List, Tuple
Restriction = Optional[List[Tuple[int, int]]]

from ._auxilliary import rotation_matrix, calcule_base
from ._exchage_map import ExchangeMap
from ._transform_molecule import find_atom_random_displ, move_mol_atom
from ._backend import (minimize_molecules, Chi2Calculator, accept_metropolis,
                       check_backend_installed)
from ._alignment import (Alignment, remove_hydrogens, guess_residue_restrains,
                         guess_protein_restrains)
from ._manager import Manager
from ._represent import (interactive_restrictions, compare_alignment,
                         compare_molecules)


__all__ = ["rotation_matrix", "calcule_base", "ExchangeMap",
           "find_atom_random_displ", "move_mol_atom", "Chi2Calculator",
           "accept_metropolis", "minimize_molecules", "Alignment",
           "interactive_restrictions", "comparate_alignment",
           "Manager", "remove_hydrogens",
           "guess_residue_restrains", "guess_protein_restrains", "Restriction"]


_DATA_DIRNAME = os.path.join(os.path.dirname(__file__), 'data')
DATA_FILES_PATH = {f: os.path.join(_DATA_DIRNAME, f)
                   for f in os.listdir(_DATA_DIRNAME)}
