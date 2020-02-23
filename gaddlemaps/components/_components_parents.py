#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
This module contains parent classes for components.

'''

from abc import abstractmethod


def are_the_same_object(object_type, *instances):
    """
    Check if the input instances correspond to the same object.

    Returns True if all instances represents same object. The inputs
    should be instance of object_type.

    Parameters
    ----------
    object_type : type
        The type of the instances.
    *instances
        Variable length list of instances instance.

    Raises
    ------
    TypeError
        If one of the input instances has wrong type.
    """

    for index, instance in enumerate(instances):
        if not isinstance(instance, object_type):
            msg = ('The instance number {} is not an instance of Instance.'
                   '').format(index)
            raise TypeError(msg)
        if not index:
            continue
        if instance != instances[0]:
            return False
    return True


class GeneralAtom(object):
    """
    The parent of every atom classes.

    """

    def __eq__(self, element):
        if isinstance(element, GeneralAtom):
            condition = ((self.resname == element.resname) and
                         (self.atomname == element.atomname))
            return condition
        return False

    def __ne__(self, element):
        return not self == element

    @staticmethod
    def are_the_same_atom(*atoms):
        """
        Check if the input atoms correspond to the same atom.

        Returns True if all atoms represents same atom. The inputs atoms should
        be instance of AtomGro or AtomItp.

        Parameters
        ----------
        *atoms
            Variable length list of atoms instance.

        Raises
        ------
        TypeError
            If one of the input atoms has wrong type.
        """

        return are_the_same_object(GeneralAtom, *atoms)


class GeneralMolecule(object):
    """
    The parent of every molecule classes.

    """

    @abstractmethod
    def __len__(self):
        pass

    def __eq__(self, element):
        if len(element) == len(self):
            for at1, at2 in zip(element, self):
                if at1 != at2:
                    return False
            return True
        return False

    def __ne__(self, element):
        return not self == element

    @staticmethod
    def are_the_same_molecule(*molecules):
        """
        Check if the input molecules correspond to the same molecule.

        Returns True if all molecules represents same molecule. The inputs
        should be instance of GeneralMolecule.

        Parameters
        ----------
        *molecules
            Variable length list of molecules instance.

        Raises
        ------
        TypeError
            If one of the input molecules has wrong type.
        """

        for index, molecule in enumerate(molecules):
            if not isinstance(molecule, GeneralMolecule):
                msg = ('The molecule number {} is not an molecule of Molecule.'
                       '').format(index)
                raise TypeError(msg)
            if not index:
                continue
            for at1, at2 in zip(molecule, molecules[0]):
                if at1 != at2:
                    return False
        return True
