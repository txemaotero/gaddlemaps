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
Molecular Simulation Components
===============================

This submodule contains useful objects to workaround with components from
molecular dynamics simulations (molecules, atoms...). You can load the
information of these components just from the coordinates files (.gro) or add
information about atom types or how they are bonded loading the topologies.

To start, you can load an object with a .gro file of a simulation system with
the class SystemGro. You can iterate through this object accessing the
residues in the system which will be instances of Residue class. At the same
time, you can iterate through the atoms (AtomGro instances) of a Residue
instance.

If you want to include the information from topology files (up to now, just .itp
files with gromacs format are compatible) to the system you can
initialize an instance of the System class. In this case, the System is
formed by Molecule objects (which are combinations of Residue and
MoleculeTop objects) and the Molecule is formed by Atom objects (combination of
AtomGro and AtomTop).

"""

from typing import List, Union

from ._components_top import MoleculeTop, AtomTop

from ._residue import Residue, AtomGro
from ._components import Atom, Molecule
from ._system import System, SystemGro


__all__ = ["AtomGro", "AtomTop", "Atom", "Residue", "MoleculeTop",
           "Molecule", "SystemGro", "System", "are_connected"]

AtomsConnected = Union[List[Atom], List[AtomTop]]

def are_connected(atoms:AtomsConnected) -> bool:
    """
    Check if the input atoms are connected.

    Parameters
    ----------
    atoms: List of AtomTop or Atom
        The list of atoms to compute if they are connected.

    Returns
    -------
    connected : bool
       True if the atoms are connected, else False.
    """

    connected_atoms: List[int] = []
    _find_connected_atoms(atoms, 0, connected_atoms)
    return len(connected_atoms) == len(atoms)


def _find_connected_atoms(atoms: AtomsConnected, index: int,
                          connected: list):
     if index not in connected:
         connected.append(index)
     for new_index in atoms[index].bonds:
         if new_index in connected:
             continue
         _find_connected_atoms(atoms, new_index, connected)
