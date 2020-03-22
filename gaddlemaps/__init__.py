# -*- coding: utf-8 -*-
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


from ._auxilliary import rotation_matrix, calcule_base
from ._exchage_map import ExchangeMap
from ._transform_molecule import find_atom_random_displ, move_mol_atom
from ._backend import minimize_molecules, Chi2Calculator, accept_metropolis
from ._alignment import (Alignment, remove_hydrogens, guess_residue_restrains,
                         guess_protein_restrains)
from ._manager import Manager


__all__ = ["rotation_matrix", "calcule_base", "ExchangeMap",
           "find_atom_random_displ", "move_mol_atom", "Chi2Calculator",
           "accept_metropolis", "minimize_molecules", "Alignment",
           "Manager", "remove_hydrogens",
           "guess_residue_restrains", "guess_protein_restrains"]
