'''
This module contains Atom and Molecule objects which connect information from
both gro and itp files.
'''

from collections import defaultdict
from typing import Any, Dict, Iterator, List, Tuple, DefaultDict

from scipy.spatial.distance import euclidean

from . import AtomGro, AtomTop, MoleculeTop, Residue


class Atom:
    """
    An atom class that wraps the AtomTop and AtomGro classes.

    You can access to the methods and attributes that both AtomGro and
    AtomItp have. To create the atom object, both input atoms should have the
    same resname and name attributes. On the other hand, only the attributes
    from the AtomGro can be changed (e.g. positions, velocities, ...) excluding
    the resname and name. 

    Parameters
    ----------
    atom_top : AtomTop
        The AtomTop object.
    atom_gro : AtomGro
        The AtomGro object.

    Raises
    ------
    IOError
        If the atom_gro and atom_top do not correspond to the same atom.
    TypeError
        If the inputs are not instances of the corresponding classes.

    """

    def __init__(self, atom_top: AtomTop, atom_gro: AtomGro):
        if not isinstance(atom_gro, AtomGro):
            raise TypeError('atom_gro input have to be an AtomGro instance.')
        if not isinstance(atom_top, AtomTop):
            raise TypeError('atom_top input have to be an AtomTop instance.')
        if atom_gro.residname != atom_top.residname:
            raise IOError((f'Input atoms do not match:\n-{atom_gro}'
                           f'\n-{atom_top}'))
        self._atom_gro = atom_gro
        self._atom_top = atom_top

    @property
    def atom_gro(self) -> AtomGro:
        """
        AtomGro : The input AtomGro object
        """
        return self._atom_gro

    @property
    def atom_top(self) -> AtomTop:
        """
        AtomTop : The input AtomTop object
        """
        return self._atom_top

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._atom_gro, attr):
            return getattr(self._atom_gro, attr)
        elif hasattr(self._atom_top, attr):
            return getattr(self._atom_top, attr)
        else:
            raise AttributeError(f'Atom object has no attribute {attr}')

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_atom_top', '_atom_gro']:
            super(Atom, self).__setattr__(attr, value)
        # Special cases that are present in both atoms
        elif attr in ('resname', 'name'):
            setattr(self._atom_top, attr, value)
            setattr(self._atom_gro, attr, value)
        elif attr in super(Atom, self).__dir__():
            super(Atom, self).__setattr__(attr, value)
        elif hasattr(self._atom_top, attr):
            setattr(self._atom_top, attr, value)
        elif hasattr(self._atom_gro, attr):
            setattr(self._atom_gro, attr, value)
        else:
            super(Atom, self).__setattr__(attr, value)

    def __str__(self) -> str:
        string = (f'Atom {self.name} of molecule {self.resname}'
                  f' with residue number {self.resid}.')
        return string

    __repr__ = __str__

    def __hash__(self) -> int:
        return hash(self._atom_top)

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Atom, self).__dir__())
        _dir.update(dir(self._atom_gro))
        _dir.update(dir(self._atom_top))
        return list(_dir)

    def copy(self) -> 'Atom':
        """
        Returns a copy of self.

        Only the gro atom is copied. The atom_top remains the same.

        Returns
        -------
        new_atom : Atom
            The copied atom.

        """
        return Atom(self._atom_top, self._atom_gro.copy())

    def __eq__(self, atom: Any) -> bool:
        if isinstance(atom, Atom):
            return self.residname == atom.residname
        return False


class Molecule:
    """
    Loads a molecule combining a MoleculeTop and a list of Residue.

    This class wraps all the features of both MoleculeTop and Residues which
    conform the molecule. When an object is initialized a copy of the input
    residues is stored (to avoid undesired attribute changes).

    Parameters
    ----------
    molecule_top : MoleculeTop
        The object with the bonds information of the molecule.
    residues : List of Residue
        An list with the residues that constitute the molecule.

    Raises
    ------
    TypeError
        If the input are instances of wrong type.
    IOError
        If residues do not constitute the molecule_top.

    """

    def __init__(self, molecule_top: MoleculeTop, residues: List[Residue]):
        if not _molecule_top_and_residues_match(molecule_top, residues):
            raise IOError(('The molecule can not be initialized. '
                           'The input residues have not the same atoms that '
                           'molecule_top has.'))

        self._molecule_top = molecule_top
        self._residues: List[Residue] = []
        self._atoms: List[Atom] = []
        # Save the residue number of each atom
        self._each_atom_resid: List[int] = []

        # Initialize attributes
        index = 0
        for res_index, res in enumerate(residues):
            self._residues.append(res.copy())
            self._each_atom_resid.append(res_index)
            for atom in res:  # type: ignore
                self._atoms.append(Atom(molecule_top[index], atom))
                index += 1

    @property
    def molecule_top(self) -> MoleculeTop:
        """
        MoleculeTop : The object with the topology information of the molecule.
        """
        return self._molecule_top

    @property
    def residues(self) -> List[Residue]:
        """
        List of Residue : a list with the Residue objects that constitute the
        molecule.
        """
        return self._residues

    def __getitem__(self, index: int) -> 'Atom':
        return self._atoms[index]

    def __len__(self) -> int:
        return len(self._atoms)

    def __iter__(self) -> Iterator['Atom']:
        for atom in self._atoms:
            yield atom

    def __getattr__(self, attr: str) -> Any:
        if hasattr(self._molecule_top, attr):
            return getattr(self._molecule_top, attr)
        raise AttributeError(f'Molecule object has no attribute {attr}')

    def __setattr__(self, attr: str, value: Any):
        if attr in ['_molecule_itp', '_residues']:
            super(Molecule, self).__setattr__(attr, value)
        elif attr in super(Molecule, self).__dir__():
            super(Molecule, self).__setattr__(attr, value)
        elif hasattr(self._molecule_top, attr):
            setattr(self._molecule_top, attr, value)
        else:
            super(Molecule, self).__setattr__(attr, value)

    def __str__(self) -> str:
        return f'Molecule of {self.name}.'

    __repr__ = __str__

    def __eq__(self, molecule: Any) -> bool:
        if (isinstance(molecule, Molecule) and
            molecule.name == self.name and
            len(molecule) == len(self)):

            for at1, at2 in zip(self, molecule):
                if at1 != at2:
                    return False
            return True
        return False

    def __ne__(self, molecule: Any) -> bool:
        return not self == molecule

    def __dir__(self) -> List[str]:
        """
        For Ipython autocompletion.
        """
        _dir = set(super(Molecule, self).__dir__())
        _dir.update(dir(self._molecule_top))
        return list(_dir)

    @property
    def atoms(self) -> List['Atom']:
        """
        list of Atom: List with the copies of the atoms of the molecule.
        """
        return [atom.copy() for atom in self._atoms]

    @property
    def resnames(self) -> List[str]:
        """
        List of string: A list with the names of the residues constituting the
            molecule.

        To set this property, a list of strings with the same lenght as the
        original must be passed. This will change each residue resname. You can
        also pass just a string and this will set all the residue names to the
        same value.
        """
        return [res.resname for res in self._residues]

    @resnames.setter
    def resnames(self, new_resnames: List[str]):
        if isinstance(new_resnames, list) and isinstance(new_resnames[0], str):
            if len(new_resnames) != len(self.resnames):
                raise ValueError(('You should provide a list with'
                                  f' {len(self.resnames)} residues instead of'
                                  f' {len(new_resnames)}.'))

            for atom, res_index in zip(self, self._each_atom_resid):
                atom.resname = new_resnames[res_index]
            for res, resname in zip(self._residues, new_resnames):
                res.resname = resname
        elif isinstance(new_resnames, str):
            for atom in self:
                atom.resname = new_resnames
            for res in self._residues:
                res.resname = new_resnames
        else:
            raise TypeError('Resnames must be a string or a list of string.')

    def copy(self, new_residues: List[Residue] = None) -> 'Molecule':
        """
        Returns a copy of the molecule.

        If new_molecule_gro is passed, the old molecule_gro will be replaced
        to update the positions. This is used in the extrapolation step.

        Parameters
        ----------
        new_residues
            List of residues to replace the original positions.

        Returns
        -------
        molecule : Molecule
            The copy of the molecule.

        """
        # Maybe this should be optimized
        if new_residues is None:
            new_residues = self._residues
        return Molecule(self._molecule_top, new_residues)

    def index(self, atom: 'Atom') -> int:
        """
        Returns the index of the atom in the molecule.

        Parameters
        ----------
        atom : Atom
            The atom to find the index.

        Returns
        -------
        index : int
            The index of the atom in the molecule.

        """
        return self._atoms.index(atom)

    @property
    def bonds_distance(self) -> Dict[int, List[Tuple[int, float]]]:
        """
        dict of int to list of tuple(int, float): A complex data structure
            that collect the information of the bond distances. The key of the
            property corresponds to the atom index in the molecule. The value
            is a list with tuples. For each tuple, the first value corresponds
            with the index of the bonded atom and the second is the length of
            the bond. This property is used in the alignment process in gaddle
            maps.
        """
        bond_info: DefaultDict[int, List[Tuple[int, float]]] = defaultdict(list)
        for index, atom in enumerate(self):
            for index_to in atom.bonds:
                atom_to = self[index_to]
                distance = euclidean(atom.position, atom_to.position)
                bond_info[index].append((index_to, distance))
        return dict(bond_info)

    @classmethod
    def from_files(cls, fgro: str, ftop: str) -> 'Molecule':
        """
        Loads the molecule from gro and a compatible topology file.

        Parameters
        ----------
        fgro : str
            The file name with the gro.
        ftop : str
            The file name with the top.

        Returns
        -------
        molecule : Molecule
            The initialized molecule.

        """
        from . import System
        syst = System(fgro, ftop)
        if len(syst) == 1:
            return syst[0]
        raise IOError('The input gro file has more than one molecule.')


def _molecule_top_and_residues_match(molecule_top: MoleculeTop,
                                     residues: List[Residue]) -> bool:
    if len(molecule_top) != sum(len(res) for res in residues):
        return False
    index = 0
    for res in residues:
        for atom in res:  # type: ignore
            if atom.residname != molecule_top[index]:
                return False
            index += 1
    return True
