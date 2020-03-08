"""
Molecular Simulation Components
===============================

This submodule contains useful objects to workaround with components from
molecular dynamics simulations (molecules, atoms...). You can load the
information of these components just from the coordinates files (.gro) or add
information  about atom types or how they  are bonded loading the topologies
(.itp).

To start, you can load an object with a .gro file of a simulation system with
the class SystemGro. You can iterate through this object accessing the molecules
in the system (fragments from the system with different residues) which will be
instances of MoleculeGro class. At the same time, you can iterate through the
atoms (AtomGro instances) of a MoleculeGro instance.

If you want to include the information from .itp files to the system you can
initialize an instance of the System class. In this case, the System is
formed by Molecule objects (which is the combination of MoleculeGro and
MoleculeItp objects) and the Molecule is formed by Atom objects (combination of
AtomGro and AtomItp).

There is also an special class (MacromoleculeGro) to build Molecule objects with
multiple residues (like proteins) as MoleculeGro is restricted to work just with
one residue.
"""

from ._components_parents import GeneralAtom, GeneralMolecule
from ._components_itp import MoleculeItp, AtomItp

from ._components_top import MoleculeTop, AtomTop

# from ._components_gro import MoleculeGro, AtomGro, MacromoleculeGro
from ._residue import Residue, AtomGro
from ._components import Atom, Molecule, InfoDict
from ._system import System, SystemGro
from ..parsers import GroLine


__all__ = ["GeneralAtom", "GeneralMolecule", "MoleculeItp", "AtomItp",
           "MoleculeGro", "AtomGro", "MacromoleculeGro", "Atom", "Molecule",
           "System", "SystemGro", "GroLine", "InfoDict"]
