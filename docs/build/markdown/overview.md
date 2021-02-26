# GADDLE Maps simplified documentation

## Basic Classes


### class gaddlemaps.Manager(system)
Class to manage the mapping process of a simulation system.

This class has methods to allow you to change the resolution of a
simulation from the .gro file of the system, the .itps of the molecules you
want to map (avoid the mapping of solvent molecules as you can resolvate the
system once it is mapped), and one .gro and one .itp for the molecules to
map in the final resolution.

This class has to be initialized with a System object but it can also be
initialized with the simulation files through self.from_files method. Then
you should specified the molecules in the final resolution using the
“add_end_molecules” method. This method will attach the molecules in the
final resolution with the corresponding one in the initial resolution
looking to pair of molecules that have the same name. If your molecules have
different names you can attach them manually by accessing the
molecule_correspondence attribute and setting the end attribute of the
Alignemnt object associated to the molecule you want to map. For example,
say you started a Manager with a system with POPC molecules and you want to
replace them with other molecules (VTE):

> Example:
> >>> vet_molecule = Molecule(vte_gro, vte_itp)
> >>> manager = Manager(System)
> >>> manager.molecule_correspondence[‘POPC’].end = vte_molecule

Once you have set the molecules in the final resolution you can call the
“align_molecules” method toe find the optimum overlap between molecules in
both resolution. Then you have to calculate the exchange maps that will be
used to extrapolate the found overlap to the rest of molecular configuration
in the system. This can be done calling the “calculate_exchange_maps”
method. Finally, you can call the “extrapolate_system” method to write a
.gro file with the system but now with the molecules in the desired final
resolution.


* **Parameters**

    **system** (*System*) – The simulation system to be mapped. It has to be an instance of System
    from system_components module.



#### molecule_correspondence()
A dictionary with the name of the loaded molecules as keys and
Alignment objects as value.


* **Type**

    dict of str: Alignment



* **Attributes**

    `complete_correspondence`

        dict of str: Alignment


### Methods

| `add_end_molecule`(molecule)

 | Add a new molecule in the end resolution to the correct Alignment.

 |
| `add_end_molecules`(\*molecules)

           | Add multiple molecules at once.

                                                                                                                                                             |
| `align_molecules`([restrictions, …])

      | Starts the alignment engine to find the optimal overlap between molecules

                                                                                                                   |
| `calculate_exchange_maps`([scale_factor])

 | Runs the alignment engine and calculate the exchange maps.

                                                                                                                                  |
| `extrapolate_system`(fgro_out)

            | Loops over the molecules in self.system and applies the exchange map.

                                                                                                                       |
| `from_files`(f_system_gro, \*ftops)

        | Build the object using the system .gro file and molecules topologies.

                                                                                                                       |
| `interactive_restrictions`([style])

       | Creates the widget to generate the restrictions of all the species in the alignment.

                                                                                                        |
| `parse_restrictions`([restrictions, …])

   | Checks the format and validates of the restrictions for the alignment.

                                                                                                                      |
<!-- !! processed by numpydoc !! -->

#### add_end_molecule(molecule)
Add a new molecule in the end resolution to the correct Alignment.


* **Parameters**

    **molecule** (*Molecule*) – The molecule in the end resolution.



* **Raises**

    
    * [**KeyError**](https://docs.python.org/3/library/exceptions.html#KeyError) – If the molecule is not found in the system.


    * [**TypeError**](https://docs.python.org/3/library/exceptions.html#TypeError) – If the molecule is not instance of Molecule


    * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If the molecule does not match with the start.


<!-- !! processed by numpydoc !! -->

#### add_end_molecules(\*molecules)
Add multiple molecules at once.


* **Parameters**

    **\*molecule** (*Molecule*) – The molecules in the end resolution.



* **Raises**

    
    * [**KeyError**](https://docs.python.org/3/library/exceptions.html#KeyError) – If the molecule is not found in the system.


    * [**TypeError**](https://docs.python.org/3/library/exceptions.html#TypeError) – If the molecule is not instance of Molecule


    * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If the molecule does not match with the start.


<!-- !! processed by numpydoc !! -->

#### align_molecules(restrictions=None, deformation_types=None, ignore_hydrogens=None, parse_restrictions=True)
Starts the alignment engine to find the optimal overlap between molecules


* **Parameters**

    
    * **restrictions** (*dict of str: list of tuple of int**, **optional*) – A dictionary with the molecules names as keys and a list of tuples
    with pairs of atom numbers corresponding to start and end atoms
    molecules. The align will be performed privileging configurations
    where those atoms are close. By default, restrictions will be set
    to [] for every molecule.

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT:** INDEX ARE REFERENCED TO THE ATOM NUMBER IN THE .itp FILE

    **NOTE**: If the restrictions were previously parsed, the
    parse_restrictions option must be set to False to avoid to parse
    the restrictions two times.



    * **deformation_types** (*dict of str: tuple of int**, **optional*) – A dictionary with the molecules names as keys and the value
    specifies the type of the minimization. Possible options:

    > 0 : Translation
    > 1 : Rotation
    > 2 : Individual atom move

    If it is None, all the possibilities are chosen for every molecule.



    * **ignore_hydrogens** (*dict of str: bool**, **optional*) – A dictionary with the molecules names as keys and a bool as value.
    If True, hydrogen atoms will not be included in the minimization
    of the distances. This will be only applied to the molecule which
    is not moved in the alignment engine. If it is None, it will be set
    as True for every molecule.


    * **parallel** (*Bool**, **optional*) – Not implemented. In the future will run the alignment for each
    molecule in parallel.


<!-- !! processed by numpydoc !! -->

#### calculate_exchange_maps(scale_factor=0.5)
Runs the alignment engine and calculate the exchange maps.


* **Parameters**

    **scale_factor** ([*float*](https://docs.python.org/3/library/functions.html#float)*, **optional*) – The compression factor to apply to mapped molecules.


<!-- !! processed by numpydoc !! -->

#### property complete_correspondence()
Alignment
A dictionary with the name of the loaded molecules as keys and
Alignment objects as value if it has start and end init.


* **Type**

    dict of str


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Dict`](https://docs.python.org/3/library/typing.html#typing.Dict)[[`str`](https://docs.python.org/3/library/stdtypes.html#str), `Alignment`]



#### extrapolate_system(fgro_out)
Loops over the molecules in self.system and applies the exchange map.


* **Parameters**

    **fgro_out** (*string*) – Gro file name to save the system in the final resolution.



* **Raises**

    [**SystemError**](https://docs.python.org/3/library/exceptions.html#SystemError) – If the exchange maps are not initialized.


<!-- !! processed by numpydoc !! -->

#### classmethod from_files(f_system_gro, \*ftops)
Build the object using the system .gro file and molecules topologies.


* **Parameters**

    
    * **f_system_gro** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Gromacs file path with the system information.


    * **\*ftops** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)*, **optional*) – A list with the topology files path of the molecules to load.



* **Returns**

    **manager** – The built mapping manager.



* **Return type**

    Manager


<!-- !! processed by numpydoc !! -->

#### interactive_restrictions(style=None)
Creates the widget to generate the restrictions of all the species in the
alignment. It generates the final representation for the widget.


* **Parameters**

    **style** (*Optional**[*[*int*](https://docs.python.org/3/library/functions.html#int)*]*) – An integer that determine which style will be used to represent the
    widget for each specie.

    > 0: One tab per specie.
    > 1: Accordion, when one specie opens the other collapse
    > 2: Vertically aligned, one over the other

    The default value is 2. This is the only one fully operational, in the
    other ones it is necessary to manually refresh the widget in the
    notebook when changing between species.




* **Return type**

    [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`ForwardRef`](https://docs.python.org/3/library/typing.html#typing.ForwardRef), [`Dict`](https://docs.python.org/3/library/typing.html#typing.Dict)[[`str`](https://docs.python.org/3/library/stdtypes.html#str), [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`int`](https://docs.python.org/3/library/functions.html#int)]]]]]



* **Returns**

    
    * **restriction_widget** (*ipywidgets.Widget*) – The widget that contains the constraint generator for all the species


    * **restrictions** (*Dict[str, List[Tuple[int, int]]]*) – The dictionary with the restrictions that will be generated by the
    widget for each specie. This will be initially empty and it will be
    filled as the widget is used.



<!-- !! processed by numpydoc !! -->

#### parse_restrictions(restrictions=None, guess_proteins=False)
Checks the format and validates of the restrictions for the alignment.


* **Parameters**

    
    * **restrictions** (*dict of str: list of tuple of int.*) – A dictionary with the molecules names as keys and a list of tuples
    with pairs of atom numbers corresponding to start and end atoms
    molecules. The align will be performed privileging configurations
    where those atoms are close. By default, restrictions will be set
    to [] for every molecule.

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT:** INDEX ARE REFERENCED TO THE ATOM NUMBER IN THE .itp FILE
    (IT USUALLY STARTS IN 1).



    * **guess_proteins** ([*bool*](https://docs.python.org/3/library/functions.html#bool)*, **optional*) – If True, restriction for proteins with more than 3 residues will be
    guessed using “guess_protein_restrain” function. This will
    overwrite the input restrains. Default False.



* **Returns**

    **new_restrictions** – The validated restrictions.

    **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
    (STARTS IN 0).**




* **Return type**

    dict of str: list of tuple of int.



* **Raises**

    
    * [**KeyError**](https://docs.python.org/3/library/exceptions.html#KeyError) – If a molecule name is not in the system.


    * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If the format is wrong or the index are not in the molecules.


<!-- !! processed by numpydoc !! -->

### class gaddlemaps.components.System(fgro, \*ftops)
Class to manage simulation systems.

A System object is formed by Molecule objects. Only the molecules
corresponding to the input ftops will be loaded.


* **Parameters**

    
    * **fgro** (*string*) – Gromacs file with the system information.


    * **\*ftops** (*string*) – Paths with the files with the bonds information to load molecules.



* **Raises**

    [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If one of the topology files do not match with any molecule in the
    system.



* **Attributes**

    `composition`

        Counter of str: int : For each molecule name (key), how many molecules there are (value).

    **fgro**


### Methods

| `add_ftop`(ftop)

                          | Adds and identifies the molecule from the ftop to the system.

                                                                                                                               |
| `add_molecule_top`(mol_top)

               | Adds a molecule to the system and find it in the gro file.

                                                                                                                                  |
<!-- !! processed by numpydoc !! -->

#### add_ftop(ftop)
Adds and identifies the molecule from the ftop to the system.

<!-- !! processed by numpydoc !! -->

#### add_molecule_top(mol_top)
Adds a molecule to the system and find it in the gro file.

<!-- !! processed by numpydoc !! -->

#### property composition()
For each molecule name (key), how many
molecules there are (value).


* **Type**

    Counter of str



* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Counter`](https://docs.python.org/3/library/typing.html#typing.Counter)[[`str`](https://docs.python.org/3/library/stdtypes.html#str)]



#### property fgro()
<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



### class gaddlemaps.components.Molecule(molecule_top, residues)
Loads a molecule combining a MoleculeTop and a list of Residue.

This class wraps all the features of both MoleculeTop and Residues which
conform the molecule. When an object is initialized a copy of the input
residues is stored (to avoid undesired attribute changes). This class
inherits from Residue so they have the same methods and properties
(although most of them are reimplemented).


* **Parameters**

    
    * **molecule_top** (*MoleculeTop*) – The object with the bonds information of the molecule.


    * **residues** (*List of Residue*) – An list with the residues that constitute the molecule.



* **Raises**

    
    * [**TypeError**](https://docs.python.org/3/library/exceptions.html#TypeError) – If the input are instances of wrong type.


    * [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If residues do not constitute the molecule_top.



* **Attributes**

    `atoms`

        list of Atom: List with the copies of the atoms of the molecule.

    `atoms_ids`

        list of int: A list with the ids of the atoms in the residues.

    `atoms_positions`

        numpy.ndarray((N, 3)) : An array with the atoms positions.

    `atoms_velocities`

        numpy.ndarray((N, 3)) or None : An array with the atoms velocities.

    `bonds_distance`

        dict of int to list of tuple(int, float): A complex data structure

    `distance_to_zero`

        float : The distance between the geometric_center and (0, 0, 0)

    `geometric_center`

        numpy.ndarray(3): Coordinates of the geometric center of the residue.

    `molecule_top`

        MoleculeTop : The object with the topology information of the molecule.

    `resid`

        int: Residue number of the residue.

    `residname`

        string: An identifier of the residue (resid+name)

    `resids`

        List of int: A list with the residues ids of the residues constituting the molecule.

    `residues`

        List of Residue : a list with the Residue objects that constitute the molecule.

    `resname`

        string: Resname of the residue.

    `resnames`

        List of string: A list with the names of the residues constituting the molecule.

    `x`

        float: The x coordinate of the geometric center of the residue.

    `y`

        float: The y coordinate of the geometric center of the residue.

    `z`

        float: The z coordinate of the geometric center of the residue.


### Methods

| `copy`([new_residues])

                    | Returns a copy of the molecule.

                                                                                                                                                             |
| `deep_copy`([new_residues])

               | Returns a deep copy of the molecule.

                                                                                                                                                        |
| `distance_to`(residue[, box_vects, inv])

  | Returns the distance between self and residue.

                                                                                                                                              |
| `from_files`(fgro, ftop)

                  | Loads the molecule from gro and a compatible topology file.

                                                                                                                                 |
| `index`(atom)

                             | Returns the index of the atom in the molecule.

                                                                                                                                              |
| `move`(displacement)

                      | Moves the residue a given displacement vector.

                                                                                                                                              |
| `move_to`(new_position)

                   | Moves the residue geometric_center to new_position.

                                                                                                                                         |
| `remove_atom`(atom)

                       | Removes a given atom from the residue.

                                                                                                                                                      |
| `rotate`(rotation_matrix)

                 | Rotate the residue around its center of mass with a given rotation matrix.

                                                                                                                  |
| `update_from_molecule_top`(mtop)

          | Modifies the Residue atoms name to match the mtop.

                                                                                                                                          |
| `write_gro`(fout)

                         | Writes a .gro file with the residue conformation.

                                                                                                                                           |
<!-- !! processed by numpydoc !! -->

#### property atoms()
List with the copies of the atoms of the molecule.


* **Type**

    list of Atom


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[`Atom`]



#### property atoms_ids()
A list with the ids of the atoms in the residues.


* **Type**

    list of int


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`int`](https://docs.python.org/3/library/functions.html#int)]



#### property atoms_positions()
An array with the atoms positions.


* **Type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)((N, 3))


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### property atoms_velocities()
An array with the atoms velocities.
If one of the atoms has no velocity this returns None.


* **Type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)((N, 3)) or [None](https://docs.python.org/3/library/constants.html#None)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)]



#### property bonds_distance()
A complex data structure
that collect the information of the bond distances. The key of the
property corresponds to the atom index in the molecule. The value
is a list with tuples. For each tuple, the first value corresponds
with the index of the bonded atom and the second is the length of
the bond. This property is used in the alignment process in gaddle
maps.


* **Type**

    dict of int to list of tuple([int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float))


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Dict`](https://docs.python.org/3/library/typing.html#typing.Dict)[[`int`](https://docs.python.org/3/library/functions.html#int), [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`float`](https://docs.python.org/3/library/functions.html#float)]]]



#### copy(new_residues=None)
Returns a copy of the molecule.

If new_molecule_gro is passed, the old residues will be replaced
to update the positions. This is used in the extrapolation step.

NOTE: With this method, the molecule_top used for the Molecule
initialization remains the same. This means that future changes in
copied molecules may affect other parts of you code. If you want a
completely independent new molecule use “deep_copy” method.


* **Parameters**

    **new_residues** ([`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`List`](https://docs.python.org/3/library/typing.html#typing.List)[`Residue`]]) – List of residues to replace the original positions.



* **Returns**

    **molecule** – The copy of the molecule.



* **Return type**

    Molecule


<!-- !! processed by numpydoc !! -->

#### deep_copy(new_residues=None)
Returns a deep copy of the molecule.

If new_molecule_gro is passed, the old residues will be replaced
to update the positions. This is used in the extrapolation step. This
method generates a new molecule that is not linked to any attribute of
the original one.


* **Parameters**

    **new_residues** ([`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`List`](https://docs.python.org/3/library/typing.html#typing.List)[`Residue`]]) – List of residues to replace the original positions.



* **Returns**

    **molecule** – The deep copy of the molecule.



* **Return type**

    Molecule


<!-- !! processed by numpydoc !! -->

#### distance_to(residue, box_vects=None, inv=False)
Returns the distance between self and residue.

residue can be a Residue instance or a 3D vector.


* **Parameters**

    
    * **residue** (*Residue** or *[*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – The residue or a point to compute the distance.


    * **box_vects** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – The box vectors to apply periodic boundary conditions.


    * **inv** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If it is True, box_vects are considered as the inverse matrix of
    the actual box_vects for a better performance.



* **Returns**

    **distance** – The euclidean distance.



* **Return type**

    [float](https://docs.python.org/3/library/functions.html#float)


<!-- !! processed by numpydoc !! -->

#### property distance_to_zero()
The distance between the geometric_center and (0, 0, 0)
without applying periodic boundary conditions.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`float`](https://docs.python.org/3/library/functions.html#float)



#### classmethod from_files(fgro, ftop)
Loads the molecule from gro and a compatible topology file.


* **Parameters**

    
    * **fgro** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The file name with the gro.


    * **ftop** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The file name with the top.



* **Returns**

    **molecule** – The initialized molecule.



* **Return type**

    Molecule


<!-- !! processed by numpydoc !! -->

#### property geometric_center()
Coordinates of the geometric center of the residue.


* **Type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)(3)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### index(atom)
Returns the index of the atom in the molecule.


* **Parameters**

    **atom** (*Atom*) – The atom to find the index.



* **Returns**

    **index** – The index of the atom in the molecule.



* **Return type**

    [int](https://docs.python.org/3/library/functions.html#int)


<!-- !! processed by numpydoc !! -->

#### property molecule_top()
The object with the topology information of the molecule.


* **Type**

    MoleculeTop


<!-- !! processed by numpydoc !! -->

* **Return type**

    `MoleculeTop`



#### move(displacement)
Moves the residue a given displacement vector.


* **Parameters**

    **displacement** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**3**)*) – An array with the displacement vector.


<!-- !! processed by numpydoc !! -->

#### move_to(new_position)
Moves the residue geometric_center to new_position.


* **Parameters**

    **new_position** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**3**)*) – An array with the new position coordinates.


<!-- !! processed by numpydoc !! -->

#### remove_atom(atom)
Removes a given atom from the residue.


* **Parameters**

    **atom** (*AtomGro*) – The atom you want to remove from the residue.


<!-- !! processed by numpydoc !! -->

#### property resid()
Residue number of the residue.


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`int`](https://docs.python.org/3/library/functions.html#int)



#### property residname()
An identifier of the residue (resid+name)


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### property resids()
A list with the residues ids of the residues

    constituting the molecule.

To set this property, a list of int with the same length as the
original must be passed. This will change each residue id. You can
also pass just an int and this will set all the residue ids to the
same value.


* **Type**

    List of int


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`int`](https://docs.python.org/3/library/functions.html#int)]



#### property residues()
a list with the Residue objects that constitute the
molecule.


* **Type**

    List of Residue


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[`Residue`]



#### property resname()
Resname of the residue.


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### property resnames()
A list with the names of the residues constituting the

    molecule.

To set this property, a list of strings with the same length as the
original must be passed. This will change each residue name. You can
also pass just a string and this will set all the residue names to the
same value.


* **Type**

    List of string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`str`](https://docs.python.org/3/library/stdtypes.html#str)]



#### rotate(rotation_matrix)
Rotate the residue around its center of mass with a given rotation
matrix.


* **Parameters**

    **rotation_matrix** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – 3x3 array with the rotation matrix.
    ADVICE: create the matrix with the help of “rotation_matrix”
    function placed in the root of the package.


<!-- !! processed by numpydoc !! -->

#### update_from_molecule_top(mtop)
Modifies the Residue atoms name to match the mtop.

This method is very useful when you have a miss-match between the
atom names in the topology and gro files. This will modify the Residue
atoms names to match the names in the topology. Make sure that all the
atoms in the topology are in the Residue.


* **Parameters**

    **mtop** (*MoleculeTop*) – The molecule to match.



* **Raises**

    [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If number of atoms in self and in mtop does not match.


<!-- !! processed by numpydoc !! -->

#### write_gro(fout)
Writes a .gro file with the residue conformation.


* **Parameters**

    **fout** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The path with the file to write the information.


<!-- !! processed by numpydoc !! -->

#### property x()
The x coordinate of the geometric center of the residue.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`float`](https://docs.python.org/3/library/functions.html#float)



#### property y()
The y coordinate of the geometric center of the residue.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`float`](https://docs.python.org/3/library/functions.html#float)



#### property z()
The z coordinate of the geometric center of the residue.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`float`](https://docs.python.org/3/library/functions.html#float)
