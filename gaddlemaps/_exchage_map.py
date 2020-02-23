#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
This module provides the functor that perform the molecule transformation
in gaddle maps.
'''

from collections import defaultdict
import numpy as np

from scipy.spatial.distance import euclidean
from .components import Molecule
from ._auxilliary import calcule_base


class ExchangeMap(object):
    """
    ExchangeMap defines a functor. It has to be initializated
    with a reference molecule and the target molecule.

    ExchangeMap(refmolecule, targetmolecule)

    Returns
    -------
        A functor that can be called with new configurations of the
        refmolecule and returns the corresponding configuration of
        the target molecule.

    Examples
    --------
       >>> BmimCG_ref = Molecule(gro,itp)
       >>> BmimAA_target = Molecule(groAA,itpAA)

       >>> transformation = ExchangeMap(BmimCG_ref, BmimAA_target)

       >>> BmimCG_new = Molecule(gro2, itp)
       >>> BmimAA_new = transformation(BmimCG_new)

    """

    def __init__(self, refmolecule, targetmolecule, scale_factor=0.5):
        self._refmolecule = refmolecule
        self._targetmolecule = targetmolecule
        self.scale_factor = scale_factor
        self._refsystems = {}
        self._equivalences = {}
        self._target_coordinates = {}
        self._calculate_refsystems(self._refmolecule)
        self._make_map()

    def _calculate_refsystems(self, molecule):
        """
        Calculate the referencce system of all the atoms that will
        be used as reference points for the coordiante
        transformation

        """
        n_atoms = len(molecule)
        # Special cases
        if n_atoms in [1, 2]:
            pos = molecule.atoms_positions
            rand_pos = [np.random.rand(3) + pos[0] for _ in range(3-n_atoms)]
            positions = np.append(pos, rand_pos, axis=0)
            coord_syst = calcule_base(positions)
            self._refsystems[hash(molecule[0])] = coord_syst
        else:
            self._calculate_refsystems_general(molecule)

    def _calculate_refsystems_general(self, molecule):
        """
        Version of _calculate_refsystems for molecules of 3 or more atoms.

        """
        # Select the atoms that are bonded at least to 2 other atoms
        for atom in molecule:
            if len(atom.bonds) >= 2:
                # Get 2 of the atoms that are bonded to the atom
                ind1, ind2 = atom.closest_atoms()
                neighbour1, neighbour2 = (molecule.hash2atom(ind1),
                                          molecule.hash2atom(ind2))
                positions = [atom.position, neighbour1.position,
                             neighbour2.position]
                # se the 3 atoms to calculate the base
                coord_syst = calcule_base(positions)
                self._refsystems[hash(atom)] = coord_syst

    def _make_map(self):
        """
        Make the conections between the atoms of the target and the references.

        It also calculates the coordinates of target in the ref coordinates
        system.

        """
        for atom in self._targetmolecule:
            name = self._find_closest_ref(atom)
            proyection = self._proyect_point(name, atom)
            self._equivalences[hash(atom)] = name
            self._target_coordinates[hash(atom)] = proyection

    def _find_closest_ref(self, targetatom):
        """
        Returns the name of the closest reference system to the atom

        """
        ref_pos = lambda index: self._refmolecule.hash2atom(index).position
        distances = [(euclidean(targetatom.position, ref_pos(index)), index)
                     for index in self._refsystems]
        return sorted(distances)[0][1]

    def _restore_point(self, atomref, proyection):
        """
        Given an atomref and a proyection restores the coordinates of the
        proyected point.

        """
        center = self._refsystems[atomref][1]
        vectores = np.array(self._refsystems[atomref][0])
        return center + np.dot(proyection, vectores)

    def _proyect_point(self, atomref, atomtarget):
        origen = self._refsystems[atomref][1]
        proyect = atomtarget.position - origen
        proyect = np.dot(self._refsystems[atomref][0], proyect)
        # Se devuelve el array escalado para evitar problemas en la minim
        return proyect * self.scale_factor

    def _restore_molecule(self):
        new_mol = self._targetmolecule.copy()
        for atom in new_mol:
            refname = self._equivalences[hash(atom)]
            proyection = self._target_coordinates[hash(atom)]
            atom.position = self._restore_point(refname, proyection)
        return new_mol

    def __call__(self, refmolecule):
        """
        This function takes as argument a molecule like the refmolecule, but in
        ther position and returns its targetmolecule equivalent as a new
        molecule with the same res_number than the input.

        """
        if not isinstance(refmolecule, Molecule):
            raise TypeError("Argument must be a Molecule")
        if self._refmolecule != refmolecule:
            raise TypeError(("refmolecule must be:\n{}"
                             "").format(self._refmolecule))

        self._calculate_refsystems(refmolecule)
        new_mol = self._restore_molecule()
        # change the indexes
        new_mol.resid = refmolecule.resid
        return new_mol

    @property
    def equivalences(self):
        """
        dict of int to list of int : {r1_atom_index: [closest_r2_atoms_indexs]}
        """
        rev = defaultdict(list)
        for ind_r1, ind_r2 in self._equivalences.items():
            rev[ind_r2].append(ind_r1)
        return dict(rev)
