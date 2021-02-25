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
This module contains the class that is used to manage the inputs in the
mapping process, starts the alignment and extrapolate the system.
'''

from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

from . import Alignment, guess_protein_restrains
from .components import Molecule, System
from .parsers import open_coordinate_file

if TYPE_CHECKING:
    from ipywidgets import Widget

Deformations = Dict[str, Optional[Tuple[int, ...]]]
Restrictions = Dict[str, Optional[List[Tuple[int, int]]]]


class Manager:
    """
    Class to manage the mapping process of a simulation system.

    This class has methods to allow you to change the resolution of a
    simulation from the .gro file of the system, the .itps of the molecules you
    want to map (avoid the mapping of solvent molecules as you can resolvate the
    system once it is mapped), and one .gro and one .itp for the molecules to
    map in the final resolution.

    This class has to be initialized with a System object but it can also be
    initialized with the simulation files through self.from_files method. Then
    you should specified the molecules in the final resolution using the
    "add_end_molecules" method. This method will attach the molecules in the
    final resolution with the corresponding one in the initial resolution
    looking to pair of molecules that have the same name. If your molecules have
    different names you can attach them manually by accessing the
    molecule_correspondence attribute and setting the end attribute of the
    Alignemnt object associated to the molecule you want to map. For example,
    say you started a Manager with a system with POPC molecules and you want to
    replace them with other molecules (VTE):

        Example:
        >>> vet_molecule = Molecule(vte_gro, vte_itp)
        >>> manager = Manager(System)
        >>> manager.molecule_correspondence['POPC'].end = vte_molecule

    Once you have set the molecules in the final resolution you can call the
    "align_molecules" method toe find the optimum overlap between molecules in
    both resolution. Then you have to calculate the exchange maps that will be
    used to extrapolate the found overlap to the rest of molecular configuration
    in the system. This can be done calling the "calculate_exchange_maps"
    method. Finally, you can call the "extrapolate_system" method to write a
    .gro file with the system but now with the molecules in the desired final
    resolution.

    Parameters
    ----------
    system : System
        The simulation system to be mapped. It has to be an instance of System
        from system_components module.

    Attributes
    ----------
    molecule_correspondence : dict of str: Alignment
        A dictionary with the name of the loaded molecules as keys and
        Alignment objects as value.

    """

    def __init__(self, system: System):
        self.system = system
        self.molecule_correspondence: Dict[str, Alignment] = {
            mol.name: Alignment(start = mol)
            for mol in self.system.different_molecules
        }

    @classmethod
    def from_files(cls, f_system_gro: str, *ftops: str) -> 'Manager':
        """
        Build the object using the system .gro file and molecules topologies.

        Parameters
        ----------
        f_system_gro : str
            Gromacs file path with the system information.
        *ftops : str, optional
            A list with the topology files path of the molecules to load.

        Returns
        -------
        manager : Manager
            The built mapping manager.

        """
        sys = System(f_system_gro, *ftops)
        return Manager(sys)

    def extrapolate_system(self, fgro_out: str):
        """
        Loops over the molecules in self.system and applies the exchange map.

        Parameters
        ----------
        fgro_out : string
            Gro file name to save the system in the final resolution.

        Raises
        ------
        SystemError
            If the exchange maps are not initialized.

        """
        complete_correspondence = self.complete_correspondence
        # Check if there is something to map
        if not complete_correspondence:
            raise SystemError('There are not loaded molecules correspondence.')
        # Check if the exchange maps are calculated
        for align in complete_correspondence.values():
            if align.exchange_map is None:
                raise SystemError(('Before extrapolating the system, '
                                   'calculate_exchange_maps method must be '
                                   'called.'))

        with open_coordinate_file(fgro_out, 'w') as fgro:
            fgro.comment = self.system.system_gro.comment_line
            fgro.box_matrix = self.system.system_gro.box_matrix
            atom_index = 1
            for mol in self.system:
                name = mol.name
                if name not in complete_correspondence:
                    continue
                new_mol = complete_correspondence[name].exchange_map(mol)  # type: ignore
                for atom in new_mol:
                    line = atom.gro_line()
                    line[3] = atom_index
                    atom_index += 1
                    fgro.writeline(line)

    @property
    def complete_correspondence(self) -> Dict[str, Alignment]:
        """
        dict of str: Alignment
            A dictionary with the name of the loaded molecules as keys and
            Alignment objects as value if it has start and end init.
        """

        return {n: ali for n, ali in self.molecule_correspondence.items()
                if (ali.end is not None) and (ali.start is not None)}

    def add_end_molecule(self, molecule: Molecule):
        """
        Add a new molecule in the end resolution to the correct Alignment.

        Parameters
        ----------
        molecule : Molecule
            The molecule in the end resolution.

        Raises
        ------
        KeyError
            If the molecule is not found in the system.
        TypeError
            If the molecule is not instance of Molecule
        ValueError
            If the molecule does not match with the start.

        """
        if not isinstance(molecule, Molecule):
            raise TypeError(('The end attribute has to be an instance of'
                             ' Molecule, not {}.').format(type(molecule)))
        name = molecule.name
        if name not in self.molecule_correspondence:
            raise KeyError(('There is not molecules with name {} in the system.'
                            ''.format(name)))
        self.molecule_correspondence[name].end = molecule

    def add_end_molecules(self, *molecules: Molecule):
        """
        Add multiple molecules at once.

        Parameters
        ----------
        *molecule : Molecule
            The molecules in the end resolution.

        Raises
        ------
        KeyError
            If the molecule is not found in the system.
        TypeError
            If the molecule is not instance of Molecule
        ValueError
            If the molecule does not match with the start.

        """
        for mol in molecules:
            self.add_end_molecule(mol)

    def align_molecules(self, restrictions: Restrictions = None,
                        deformation_types: Deformations = None,
                        ignore_hydrogens: Dict[str, bool] = None,
                        parse_restrictions: bool = True):
        """
        Starts the alignment engine to find the optimal overlap between molecules

        Parameters
        ----------
        restrictions : dict of str: list of tuple of int, optional
            A dictionary with the molecules names as keys and a list of tuples
            with pairs of atom numbers corresponding to start and end atoms
            molecules. The align will be performed privileging configurations
            where those atoms are close. By default, restrictions will be set
            to [] for every molecule.

            Example
            -------

            >>> restrictions = [(1, 3), (4, 5)]

            **IMPORTANT:** INDEX ARE REFERENCED TO THE ATOM NUMBER IN THE .itp FILE

            .. note:: If the restrictions were previously parsed, the
                parse_restrictions option must be set to False to avoid to parse
                the restrictions two times.

        deformation_types : dict of str: tuple of int, optional
            A dictionary with the molecules names as keys and the value
            specifies the type of the minimization. Possible options:
                0 : Translation
                1 : Rotation
                2 : Individual atom move
            If it is None, all the possibilities are chosen for every molecule.
        ignore_hydrogens : dict of str: bool, optional
            A dictionary with the molecules names as keys and a bool as value.
            If True, hydrogen atoms will not be included in the minimization
            of the distances. This will be only applied to the molecule which
            is not moved in the alignment engine. If it is None, it will be set
            as True for every molecule.
        parallel : Bool, optional
            Not implemented. In the future will run the alignment for each
            molecule in parallel.

        """

        if parse_restrictions or restrictions is None:
            restrictions = self.parse_restrictions(restrictions)
        assert isinstance(restrictions, dict)  # for typing
        deformation_types = self._parse_deformations(deformation_types)
        ignore_hydrogens = self._parse_ignore_hydrogens(ignore_hydrogens)
        mols_corr = self.complete_correspondence
        for name in restrictions:
            restr = restrictions[name]
            defor = deformation_types[name]
            ignor = ignore_hydrogens[name]
            print('Aligning {}:\n'.format(name))
            mols_corr[name].align_molecules(restr, defor, ignor)

    def parse_restrictions(self, restrictions: Restrictions = None,
                           guess_proteins: bool = False):
        """
        Checks the format and validates of the restrictions for the alignment.

        Parameters
        ----------
        restrictions : dict of str: list of tuple of int.
            A dictionary with the molecules names as keys and a list of tuples
            with pairs of atom numbers corresponding to start and end atoms
            molecules. The align will be performed privileging configurations
            where those atoms are close. By default, restrictions will be set
            to [] for every molecule.

            Example
            -------
            >>> restrictions = [(1, 3), (4, 5)]

            **IMPORTANT:** INDEX ARE REFERENCED TO THE ATOM NUMBER IN THE .itp FILE
            (IT USUALLY STARTS IN 1).
        guess_proteins : bool, optional
            If True, restriction for proteins with more than 3 residues will be
            guessed using "guess_protein_restrain" function. This will
            overwrite the input restrains. Default False.

        Returns
        -------
        new_restrictions : dict of str: list of tuple of int.
            The validated restrictions.

            **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
            (STARTS IN 0).**

        Raises
        ------
        KeyError
            If a molecule name is not in the system.
        ValueError
            If the format is wrong or the index are not in the molecules.

        """
        if restrictions is None:
            new_restrictions: Dict[str, Optional[List[Tuple[int, int]]]] = {}
            for name in self.complete_correspondence:
                new_restrictions[name] = None
                if guess_proteins:
                    start = self.molecule_correspondence[name].start
                    assert isinstance(start, Molecule)
                    resnames = start.resnames
                    if len(resnames) > 3:
                        end = self.molecule_correspondence[name].end
                        assert isinstance(end , Molecule)
                        restr = guess_protein_restrains(start, end)
                        new_restrictions[name] = restr
                        continue
            return new_restrictions
        complete_correspondence = self.complete_correspondence
        for name_comp in restrictions:
            if name_comp not in complete_correspondence:
                raise KeyError('There are no molecules in the system with name'
                               ': {}'.format(name_comp))

        new_restrictions = {}
        for name in complete_correspondence:
            if guess_proteins:
                start = self.molecule_correspondence[name].start
                assert isinstance(start, Molecule)
                resnames = start.resnames
                if len(resnames) > 3:
                    end = self.molecule_correspondence[name].end
                    assert isinstance(end, Molecule)
                    restr = guess_protein_restrains(start, end)
                    new_restrictions[name] = restr
                    continue
            if name in restrictions:
                restriction = restrictions[name]
                if not restriction:
                    new_restrictions[name] = None
                else:
                    new_restrictions[name] = self._validate_index(restriction,
                                                                  name)
            else:
                new_restrictions[name] = None
        return new_restrictions

    def _validate_index(self, restriction: List[Tuple[int, int]],
                        name: str) -> List[Tuple[int, int]]:
        """Validates the input restrictions and change from atomid to index.
        """
        mol_start = self.molecule_correspondence[name].start
        assert isinstance(mol_start, Molecule)
        mol_end = self.molecule_correspondence[name].end
        assert isinstance(mol_end, Molecule)

        msg = f'{name}: '
        msg_format = msg + ('The input restriction has wrong format. See '
                            'documentation of align_molecules method.')

        msg_index = msg + ('Error accessing the {}th atom in the {}'
                           ' resolution molecule.')
        list_res = []
        for tup in restriction:
            if not hasattr(tup, '__len__') or len(tup) != 2:
                raise ValueError(msg_format)
            try:
                ind1 = mol_start[tup[0]]
            except IndexError:
                raise ValueError(msg_index.format(tup[0], 'initial'))
            try:
                ind2 = mol_end[tup[1]]
            except IndexError:
                raise ValueError(msg_index.format(tup[1], 'final'))
            list_res.append(tup)
        return list_res

    def _parse_deformations(self, deformations: Deformations = None) -> Dict[str, Any]:
        if deformations is None:
            return {name: None for name in self.complete_correspondence}
        complete_correspondence = self.complete_correspondence
        for name in deformations:
            if name not in complete_correspondence:
                raise KeyError('There are no molecules with names {} in the '
                               'system.'.format(', '.join(deformations.keys())))
        new_def: Deformations = {}
        for name in complete_correspondence:
            if name in deformations:
                deformation = deformations[name]
                if not deformation:
                    new_def[name] = None
                else:
                    if not hasattr(deformation, '__len__'):
                        raise ValueError(('Wrong format for deformation_types.'
                                          ' See documentation of '
                                          'align_molecules method.'))
                    if not 1 <= len(deformation) <= 3:
                        raise ValueError(('Wrong format for deformation_types.'
                                          ' See documentation of '
                                          'align_molecules method.'))
                    new_def[name] = deformation
            else:
                new_def[name] = None
        return new_def

    def _parse_ignore_hydrogens(self, ignore_hydrogens: Dict[str, bool] = None) -> Dict[str, bool]:
        if ignore_hydrogens is None:
            return {name: True for name in self.complete_correspondence}
        complete_correspondence = self.complete_correspondence
        for name in ignore_hydrogens:
            if name not in complete_correspondence:
                raise KeyError(('There are no molecules with names {} in the '
                                'system.'
                                '').format(', '.join(ignore_hydrogens.keys())))
        new_ign = {}
        for name in complete_correspondence:
            if name not in ignore_hydrogens:
                new_ign[name] = True
            else:
                val = ignore_hydrogens[name]
                if not isinstance(val, bool):
                    raise ValueError(('Wrong format for ignore_hydrogens. See '
                                      'documentation of align_molecules method.'
                                      ''))
                new_ign[name] = val
        return new_ign

    def calculate_exchange_maps(self, scale_factor: float = 0.5):
        """
        Runs the alignment engine and calculate the exchange maps.

        Parameters
        ----------
        scale_factor : float, optional
            The compression factor to apply to mapped molecules.

        """
        complete_correspondence = self.complete_correspondence
        for name in complete_correspondence:
            complete_correspondence[name].init_exchange_map(scale_factor)

    def interactive_restrictions(self, style:int=None) -> Tuple['Widget',
                                                                Restrictions]:
        """
        Creates the widget to generate the restrictions of all the species in the
        alignment. It generates the final representation for the widget.

        Parameters
        ----------
        style: Optional[int]
            An integer that determine which style will be used to represent the
            widget for each specie.
                0: One tab per specie.
                1: Accordion, when one specie opens the other collapse
                2: Vertically aligned, one over the other
            The default value is 2. This is the only one fully operational, in the
            other ones it is necessary to manually refresh the widget in the
            notebook when changing between species.

        Returns
        -------
        restriction_widget: ipywidgets.Widget
            The widget that contains the constraint generator for all the species
        restrictions: Dict[str, List[Tuple[int, int]]]
            The dictionary with the restrictions that will be generated by the
            widget for each specie. This will be initially empty and it will be
            filled as the widget is used.

        """
        from ._represent import interactive_restrictions

        return interactive_restrictions(self.molecule_correspondence,
                                        style=style)