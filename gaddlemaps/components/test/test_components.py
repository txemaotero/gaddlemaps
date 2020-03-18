"""
Tests for _components submodule.
"""

from typing import Callable, List

import os
import pytest
import numpy as np

from gaddlemaps.parsers import GroFile
from gaddlemaps.components import (AtomGro, AtomTop, Atom, Molecule,
                                   MoleculeTop, Residue)
from gaddlemaps.components._components import _molecule_top_and_residues_match


AtomGenerator = Callable[..., AtomGro]

ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]


@pytest.fixture
def atom_gro_generator() -> AtomGenerator:
    """
    A function to create AtomGro instances. 
    """
    def generate_atom_gro(resid: int = 1, resname: str = 'BMIM',
                          name: str = 'N1', atomid: int = 1) -> AtomGro:
        info = (resid, resname, name, atomid, 4.668,
                3.571, 8.232, -0.2489, 0.2514, 0.1046)
        return AtomGro(info)
    return generate_atom_gro

@pytest.fixture
def atom_n1_bmim(atom_gro_generator: AtomGenerator) -> Atom:
    """
    N1 Atom instance of BMIM.
    """
    atom_gro = atom_gro_generator()
    atom_top = AtomTop('N1', 'BMIM', 1, 0)
    return Atom(atom_top, atom_gro)


@pytest.fixture
def molecule_top_bmim() -> MoleculeTop:
    """
    BMIM MoleculeTop instance.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.itp')
    return MoleculeTop(fname)


@pytest.fixture
def residue_bmim() -> Residue:
    """
    BMIM Residue instance.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/BMIM_AA.gro')
    with GroFile(fname) as fgro:
        atoms = [AtomGro(line) for line in fgro]
    return Residue(atoms)


@pytest.fixture
def molecule_bmim(molecule_top_bmim: MoleculeTop,
                  residue_bmim: Residue) -> Molecule:
    """
    BMIM Molecule instance.
    """
    return Molecule(molecule_top_bmim, [residue_bmim])


@pytest.fixture
def molecule_top_protein() -> MoleculeTop:
    """
    MoleculeTop instance for molecule with multiple residues.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/protein_CG.itp')
    return MoleculeTop(fname)


@pytest.fixture
def residue_protein() -> List[Residue]:
    """
    List of Residue instances for a molecule with multiple residues.
    """
    fname = os.path.join(ACTUAL_PATH, '../../data/protein_CG.gro')
    residues: List[Residue] = []
    with GroFile(fname) as fgro:
        atom = AtomGro(next(fgro))
        residname = atom.residname
        atoms = [atom]
        for line in fgro:
            atom = AtomGro(line)
            if residname != atom.residname:
                residues.append(Residue(atoms))
                residname = atom.residname
                atoms = [atom]
            else:
                atoms.append(atom)
    residues.append(Residue(atoms))
    return residues


@pytest.fixture
def molecule_protein(molecule_top_protein: MoleculeTop,
                     residue_protein: List[Residue]) -> Molecule:
    """
    Molecule instance with multiple residues.
    """
    return Molecule(molecule_top_protein, residue_protein)


def test_molecule_top_and_residues_match(molecule_top_bmim: MoleculeTop,
                                         residue_bmim: Residue,
                                         molecule_top_protein: MoleculeTop,
                                         residue_protein: List[Residue]):
    """
    Test the match between residues and MoleculeTop.
    """
    # different length
    assert not _molecule_top_and_residues_match(molecule_top_bmim,
                                                residue_protein)
                                                
    assert _molecule_top_and_residues_match(molecule_top_bmim,
                                            [residue_bmim])
    assert _molecule_top_and_residues_match(molecule_top_protein,
                                            residue_protein)


class TestAtom:
    """
    Wraps all the tests for the Atom class.
    """
    def test_initialization(self, atom_gro_generator: AtomGenerator):
        """
        Tests Atom initialization and own attributes
        """
        atom_gro = atom_gro_generator()
        atom_top = AtomTop('N1', 'BMIM', 1, 0)
        atom = Atom(atom_top, atom_gro)
        assert atom._atom_gro is atom_gro
        assert atom._atom_top is atom_top
        new_atom = atom.copy()
        assert new_atom == atom
        assert new_atom is not atom
        assert hash(atom) == 0

        with pytest.raises(TypeError, match="atom_gro.*"):
            atom = Atom(atom_top, atom_top)  # type:ignore
        with pytest.raises(TypeError, match="atom_top.*"):
            atom = Atom(atom_gro, atom_gro)  # type:ignore

        atom_gro.name = 'B5'
        with pytest.raises(IOError):
            atom = Atom(atom_top, atom_gro)

        atom_gro.name = 'N1'
        atom_gro.resname = 'BF4'
        with pytest.raises(IOError):
            atom = Atom(atom_top, atom_gro)

    def test_get_set_attributes(self, atom_n1_bmim: Atom):
        """
        Tests get and set attributes for Atom.
        """
        # Name
        assert atom_n1_bmim.name == 'N1'
        atom_n1_bmim.name = 'B5'
        assert atom_n1_bmim.name == 'B5'
        assert atom_n1_bmim._atom_gro.name == 'B5'
        assert atom_n1_bmim._atom_top.name == 'B5'

        # Resname
        assert atom_n1_bmim.resname == 'BMIM'
        atom_n1_bmim.resname = 'BF4'
        assert atom_n1_bmim.resname == 'BF4'
        assert atom_n1_bmim._atom_gro.resname == 'BF4'
        assert atom_n1_bmim._atom_top.resname == 'BF4'

        # resids
        with pytest.raises(AttributeError):
            resid = atom_n1_bmim.resid
        with pytest.raises(AttributeError):
            atom_n1_bmim.resid = 1

        assert atom_n1_bmim.gro_resid == 1
        atom_n1_bmim.gro_resid = 3
        assert atom_n1_bmim.gro_resid == 3
        assert atom_n1_bmim._atom_gro.resid == 3

        assert atom_n1_bmim.top_resid == 1
        atom_n1_bmim.top_resid = 5
        assert atom_n1_bmim.top_resid == 5
        assert atom_n1_bmim._atom_top.resid == 5

        # Residname is the gro one
        assert atom_n1_bmim.residname == '3BF4'

        # index
        assert atom_n1_bmim.index == 0
        assert hash(atom_n1_bmim) == 0
        atom_n1_bmim.index = 5
        assert atom_n1_bmim.index == 5
        assert hash(atom_n1_bmim) == 5
        assert atom_n1_bmim._atom_top.index == 5
        
        # bonds 
        assert atom_n1_bmim.bonds == atom_n1_bmim._atom_top.bonds
        atom_n1_bmim.bonds = {1, 32, 4}
        assert atom_n1_bmim.bonds == {1, 32, 4}
        assert atom_n1_bmim.closest_atoms() == [1, 4]
        atom_n1_bmim.connect(atom_n1_bmim)
        assert atom_n1_bmim.bonds == {1, 32, 4, 5}

        # atomid
        assert atom_n1_bmim.atomid == 1
        atom_n1_bmim.atomid = 6
        assert atom_n1_bmim.atomid == 6
        assert atom_n1_bmim._atom_gro.atomid == 6

        # element
        assert atom_n1_bmim.element == 'B'

        # position
        assert np.isclose(atom_n1_bmim.position, [4.668, 3.571, 8.232]).all()
        atom_n1_bmim.position = np.array([0, 0, 0])
        assert np.isclose(atom_n1_bmim.position, [0, 0, 0]).all()
        assert np.isclose(atom_n1_bmim._atom_gro.position, [0, 0, 0]).all()

        # velocity
        assert np.isclose(atom_n1_bmim.velocity, [-0.2489, 0.2514, 0.1046]).all()
        atom_n1_bmim.velocity = None
        assert atom_n1_bmim.velocity is None
        assert atom_n1_bmim._atom_gro.velocity is None

        # gro_line
        assert atom_n1_bmim.gro_line() == atom_n1_bmim._atom_gro.gro_line()


class TestMolecule:
    """
    Wraps all the tests for the Molecule class.
    """
    def test_init_one_residue(self, molecule_top_bmim: MoleculeTop,
                              residue_bmim: Residue):
        """
        Tests Molecule initialization with only one residue and some basic
        attributes and methods.
        """
        molec = Molecule(molecule_top_bmim, [residue_bmim])

        assert len(molec) == len(molecule_top_bmim)
        assert len(molec) == len(residue_bmim)
        assert len(molec) == 25

        for atom, atom_top, atom_gro in zip(molec, molecule_top_bmim,
                                            residue_bmim):
            atom_compare = Atom(atom_top, atom_gro)
            assert atom == atom_compare

        for atom, atom_comp in zip(molec, molec.atoms):
            assert atom == atom_comp

        assert molec[-1] == atom_compare
        assert atom_compare != molec[0]
        assert molec.index(atom_compare) == 24

        assert molec._each_atom_resid == [0] * len(molec)

        assert molec == molec
        assert molec is molec

        assert molec.copy() == molec
        assert molec.copy() is not molec

        assert len(molec.residues) == 1
        assert molec.residues[0] == residue_bmim
        assert molec.residues[0] is not residue_bmim
        assert molec.molecule_top is molecule_top_bmim

    def test_init_multi_residue(self, molecule_top_protein: MoleculeTop,
                                residue_protein: List[Residue]):
        """
        Tests Molecule initialization with multiple residues and some basic
        attributes and methods.
        """
        molec = Molecule(molecule_top_protein, residue_protein)

        assert len(molec) == len(molecule_top_protein)
        assert len(molec) == sum(len(res) for res in residue_protein)
        assert len(molec) == 34

        index = 0
        for res in residue_protein:
            for atom_gro in res:
                atom_compare = Atom(molecule_top_protein[index], atom_gro)
                assert molec[index] == atom_compare
                index += 1

        for atom, atom_comp in zip(molec, molec.atoms):
            assert atom == atom_comp

        assert molec[-1] == atom_compare
        assert molec.index(atom_compare) == 33

        assert molec._each_atom_resid == [0, 0, 0, 1, 1, 1, 2, 2, 3,
                                          3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7,
                                          8, 8, 8, 9, 9, 10, 10, 10, 11, 11,
                                          11, 12, 12]

        assert molec == molec
        assert molec is molec

        assert molec.copy() == molec
        assert molec.copy() is not molec

        assert len(molec.residues) == 13
        for res1, res2 in zip(molec.residues, residue_protein):
            assert res1 == res2
            assert res1 is not res2
        assert molec.molecule_top is molecule_top_protein

    def test_init_from_files(self):
        """
        Test the molecule initialization from files.
        """
        # TODO: Implement when the tests for System are implemented.
        pass

    def test_resnames(self, molecule_protein: Molecule,
                      molecule_bmim: Molecule):
        """
        Tests for the resnames property
        """
        assert molecule_bmim.resnames == ['BMIM']
        molecule_bmim.resnames = ['BF4']
        assert molecule_bmim.resnames == ['BF4']
        assert molecule_bmim.molecule_top.resnames == ['BF4']
        assert molecule_bmim.residues[0].resname == 'BF4'
        for atom in molecule_bmim:
            assert atom.resname == 'BF4'

        molecule_bmim.resnames = 'test'
        assert molecule_bmim.resnames == ['test']
        assert molecule_bmim.molecule_top.resnames == ['test']
        assert molecule_bmim.residues[0].resname == 'test'
        for atom in molecule_bmim:
            assert atom.resname == 'test'

        resnames_compare = ['ARG', 'ARG', 'LEU', 'LYS', 'ARG', 'LEU', 'LEU',
                            'ARG', 'ARG', 'LEU', 'LYS', 'ARG', 'LEU']
        assert molecule_protein.resnames == resnames_compare
        resnames_compare[3] = 'Test'
        molecule_protein.resnames = resnames_compare
        assert molecule_protein.resnames == resnames_compare
        for res, res_test in zip(molecule_protein.residues, resnames_compare):
            assert res.resname == res_test
        for index, atom in enumerate(molecule_protein):
            res_index = molecule_protein._each_atom_resid[index]
            assert atom.resname == resnames_compare[res_index]

        molecule_protein.resnames = 'hello'
        assert molecule_protein.resnames == ['hello'] * len(resnames_compare)
        for res, res_test in zip(molecule_protein.residues, resnames_compare):
            assert res.resname == 'hello'
        for atom in molecule_protein:
            assert atom.resname == 'hello'

    def test_get_set_attributes(self, molecule_top_bmim: MoleculeTop,
                                molecule_bmim: Molecule):
        """
        Tests get and set attributes for Molecule. Molecule getattr just access
        the molecule_top attributes.
        """
        # Name
        assert molecule_top_bmim.name == molecule_bmim.name
        molecule_bmim.name = 'Test'
        assert molecule_bmim.name == 'Test'
        assert molecule_bmim.molecule_top.name == 'Test'

        # Resname len list
        assert molecule_bmim.resname_len_list == [('BMIM', 25)]

    def test_bonds_distance(self, molecule_bmim: Molecule):
        """
        Test bonds distances for bmim.
        """
        bonds_dist = molecule_bmim.bonds_distance
        for index, atom in enumerate(molecule_bmim):
            assert set(val[0] for val in bonds_dist[index]) == atom.bonds
            for index_to, dist in bonds_dist[index]:
                diff = atom.position - molecule_bmim[index_to].position
                dist_test = np.linalg.norm(diff)
                assert round(dist, 5) == round(dist_test, 5)
