#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
This module contains the class that is used to manage the inputs in the
mapping process, starts the alignment and extrapolate the maps.
'''

from .components import System, Molecule
from .parsers import GroFile
from . import Alignment, guess_protein_restrains


class GaddleMapsManager(object):
    """
    Class to manage a simulation system mapping process.

    The class has to be init with a system object but is can also be
    initialized with the simulation files through self.from_files method.

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
    def __init__(self, system):
        self.system = system
        self.molecule_correspondence = {
            mol.name: Alignment(start=mol)
            for mol in self.system.different_molecules
        }

    @classmethod
    def from_files(cls, f_system_gro, *fitps):
        """
        Build the object using the system .gro file and molecules .itp.

        Parameters
        ----------
        f_system_gro : str
            Gromacs file path with the system information.
        *fitps : str, optional
            A list with the .itp files path of the molecules to load.

        Returns
        -------
        manager : GaddleMapsManager
            The built mapping manager.

        """
        sys = System(f_system_gro, *fitps)
        return GaddleMapsManager(sys)

    def extrapolate_system(self, fgro_out):
        """
        Loops over the molecules in self.system and applies the exchange map.

        The mapped molecules are written in fgro_out file instead of crating
        new system.

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

        with GroFile(fgro_out, 'w') as fgro:
            fgro.comment = self.system.system_gro.comment_line
            fgro.box_matrix = self.system.system_gro.box_matrix
            atom_index = 1
            for mol in self.system:
                name = mol.name
                if name not in complete_correspondence:
                    continue
                new_mol = complete_correspondence[name].exchange_map(mol)
                for atom in new_mol:
                    line = atom.gro_line()
                    line[3] = atom_index
                    atom_index += 1
                    fgro.writeline(line)

    @property
    def info(self):
        """
        dict : A dictionary with the needed information to restore the
            manager.
        """
        mol_corr_info = {n: ali.info
                         for n, ali in self.molecule_correspondence.items()}
        info = {
            'system': self.system.info,
            'molecule_correspondence': mol_corr_info,
        }

        return info
        raise AttributeError(('GaddleMapsManager object has not '
                              'info property.'))

    @classmethod
    def from_info(cls, info_dict):
        """
        Builds the manager from the information returned by info property.

        Parameters
        ----------
        info_dict : A dictionary with the needed information to restore the
            manager.

        Returns
        -------
        manager : GaddleMapsManager
            The loaded system.
        """
        man = GaddleMapsManager(System.from_info(info_dict['system']))
        mol_corr = {n: Alignment.from_info(ali)
                    for n, ali in info_dict['molecule_correspondence'].items()}
        man.molecule_correspondence = mol_corr
        return man

    @property
    def complete_correspondence(self):
        """
        dict of str: Alignment
            A dictionary with the name of the loaded molecules as keys and
            Alignment objects as value if it has start and end init.
        """

        return {n: ali for n, ali in self.molecule_correspondence.items()
                if (ali.end is not None) and (ali.start is not None)}

    def add_end_molecule(self, molecule):
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
            raise TypeError(('The start attribute has to be an instance of'
                             ' Molecule, not {}.').format(type(molecule)))
        name = molecule.name
        if name not in self.molecule_correspondence:
            raise KeyError(('There is not molecules with name {} in the system.'
                            ''.format(name)))
        self.molecule_correspondence[name].end = molecule

    def add_end_molecules(self, *molecules):
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

    def align_molecules(self, restrictions=None, deformation_types=None,
                        ignore_hydrogens=None, parallel=False,
                        parse_restrictions=True):
        """
        Starts the alignment engine to find the optimal overlap between molecs.

        Parameters
        ----------
        restrictions : dict of str: list of tuple of int, optional
            A dictionary with the molecules names as keys and a list of tuples
            with pairs of atom numbers corresponding to start and end atoms
            molecules. The align will be performed privileging configurations
            where those atoms are close. By default, restrictions will be set
            to [] for every molecule.

            Example:
            >>> restrictions = [(1, 3), (4, 5)]

            IMPORTANT: INDEX ARE REFERENCED TO THE ATOM NUMBER IN THE .itp FILE

            NOTE: If the restrictions were previously parsed, the
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

        Raises
        ------
        NotImplementedError
            If parallel is set to True.

        """

        if parallel:
            raise NotImplementedError('The parallelization is not implemented.')
        if parse_restrictions:
            restrictions = self.parse_restrictions(restrictions)
        deformation_types = self._parse_deformations(deformation_types)
        ignore_hydrogens = self._parse_ignore_hydrogens(ignore_hydrogens)
        mols_corr = self.complete_correspondence
        for name in restrictions:
            restr = restrictions[name]
            defor = deformation_types[name]
            ignor = ignore_hydrogens[name]
            print('Aligning {}:\n'.format(name))
            mols_corr[name].align_molecules(restr, defor, ignor)

    def parse_restrictions(self, restrictions, guess_proteins=False):
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

            Example:
            >>> restrictions = [(1, 3), (4, 5)]

            IMPORTANT: INDEX ARE REFERENCED TO THE ATOM NUMBER IN THE .itp FILE
            (STARTS IN 1).
        guess_proteins : bool, optional
            If True, restriction for proteins with more than 3 residues will be
            guessed using "guess_protein_restrain" function. This will
            overwrite the input restrains. Default False.

        Returns
        -------
        new_restrictions : dict of str: list of tuple of int.
            The validated restrictions.

            IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
            (STARTS IN 0).

        Raises
        ------
        KeyError
            If a molecule name is not in the system.
        ValueError
            If the format is wrong or the index are not in the molecules.

        """
        if restrictions is None:
            new_restrictions = {}
            for name in self.complete_correspondence:
                new_restrictions[name] = None
                if guess_proteins:
                    start = self.molecule_correspondence[name].start
                    resnames = start.resnames
                    if len(resnames) > 3:
                        end = self.molecule_correspondence[name].end
                        restr = guess_protein_restrains(start, end)
                        new_restrictions[name] = restr
                        continue
            return new_restrictions
        complete_correspondence = self.complete_correspondence
        for name in restrictions:
            if name not in complete_correspondence:
                raise KeyError('There are no molecules in the system with name'
                               ': {}'.format(name))

        new_restrictions = {}
        for name in complete_correspondence:
            if guess_proteins:
                start = self.molecule_correspondence[name].start
                resnames = start.resnames
                if len(resnames) > 3:
                    end = self.molecule_correspondence[name].end
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

    def _validate_index(self, restriction, name):
        """Validates the input restrictions and change from atomid to index.
        """
        mol_start = self.molecule_correspondence[name].start
        mol_end = self.molecule_correspondence[name].end

        msg_format = ('{}: The input restriction has wrong format. See '
                      'documentation of align_molecules method.')

        list_res = []
        for tup in restriction:
            if not hasattr(tup, '__len__'):
                raise ValueError(msg_format.format(name))
            if len(tup) != 2:
                raise ValueError(msg_format.format(name))
            try:
                ind1 = mol_start.hash2index(tup[0])
            except KeyError:
                msg = '{}: '.format(name)
                msg += ('Molecule in initial resolution has no atom'
                        ' with number {}'.format(tup[0]))
                raise ValueError(msg)
            try:
                ind2 = mol_end.hash2index(tup[1])
            except KeyError:
                msg = '{}: '.format(name)
                msg += ('Molecule in final resolution has no atom'
                        ' with number {}'.format(tup[1]))
                raise ValueError(msg)
            list_res.append((ind1, ind2))
        return list_res

    def _parse_deformations(self, deformations):
        if deformations is None:
            return {name: None for name in self.complete_correspondence}
        complete_correspondence = self.complete_correspondence
        for name in deformations:
            if name not in complete_correspondence:
                raise KeyError('There are no molecules with names {} in the '
                               'system.'.format(', '.join(deformations.keys())))
        new_def = {}
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

    def _parse_ignore_hydrogens(self, ignore_hydrogens):
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

    def calculate_exchange_maps(self, scale_factor=0.5):
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
