# gaddlemaps.components package

## Module contents

### Molecular Simulation Components

This submodule contains useful objects to workaround with components from
molecular dynamics simulations (molecules, atoms…). You can load the
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

<!-- !! processed by numpydoc !! -->

### class gaddlemaps.components.AtomGro(parsed_gro_line)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

It contains the information of a .gro line corresponding to an atom.

You can add atoms with the same residue name and number to form a Residue
instance.


* **Parameters**

    **parsed_gro_line** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – A list returned by the GroFile.read_line method when a correct
    formatted .gro line is input.



#### resid()
Residue number of the atom.


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)



#### resname()
Residue name of the atom.


* **Type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)



#### name()
Name of the atom.


* **Type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)



#### atomid()
Atom index in the .gro file.


* **Type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)



#### position()
The coordinates of the atom.


* **Type**

    np.ndarray(3)



#### velocity()
The velocity of the atom if it is specified in the .gro line. If not
this attribute is set to None.


* **Type**

    np.ndarray(3) or [None](https://docs.python.org/3/library/constants.html#None)



* **Attributes**

    `element`

        string : The element of the atom. It is obtained removing the non

    `residname`

        string : An identifier that contains the residue name and number. This


### Methods

| `copy`(self)

 | Returns a copy of the atom.

 |
| `gro_line`(self, parsed)

                      | Returns the gro line corresponding to the atom.

                                                                                                                                             |
<!-- !! processed by numpydoc !! -->

#### copy(self)
Returns a copy of the atom.

You can safely change the attributes of the returned atom without
changing the original one.


* **Returns**

    **new_atom** – The copied atom.



* **Return type**

    AtomGro


<!-- !! processed by numpydoc !! -->

#### property element()
The element of the atom. It is obtained removing the non
alphabetic character in the atom name.


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### gro_line(self, parsed: bool = True)
Returns the gro line corresponding to the atom.


* **Parameters**

    **parsed** ([*bool*](https://docs.python.org/3/library/functions.html#bool)*, **optional*) – If True (default), the line is returned as FileGro.read_line
    method output. Else, line with the correct .gro format is
    returned.



* **Returns**

    **gro_line** – The corresponding .gro line.



* **Return type**

    List of ([str](https://docs.python.org/3/library/stdtypes.html#str), [int](https://docs.python.org/3/library/functions.html#int) or [float](https://docs.python.org/3/library/functions.html#float)) or [str](https://docs.python.org/3/library/stdtypes.html#str)


<!-- !! processed by numpydoc !! -->

#### property residname()
An identifier that contains the residue name and number. This
should be unique for ear Residue object in a simulation system.


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



### class gaddlemaps.components.AtomTop(name, resname, resid, index)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Atom with information of its name, residue name and bonds.

It is also needed the index of the atom in the molecule. The bonds are
initialized as empty set.


* **Parameters**

    
    * **name** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Atom name


    * **resname** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Residue name of the atom.


    * **resid** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Residue index of the atom.


    * **index** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Atom index in the molecule.



#### bonds()
A set with the hash of the atoms that are connected to self.


* **Type**

    set of int



* **Attributes**

    `residname`

        string: An identifier of the residue (resid+name)


### Methods

| `closest_atoms`(self, natoms)

                 | Returns a list with natoms index of bonded atoms to self.

                                                                                                                                   |
| `connect`(self, atom)

                         | Connects self with other atom setting the bond.

                                                                                                                                             |
| `copy`(self)

                                  | Returns a copy of the current atom.

                                                                                                                                                         |
<!-- !! processed by numpydoc !! -->

#### closest_atoms(self, natoms: int = 2)
Returns a list with natoms index of bonded atoms to self.

If more than natoms are bonded self, the natoms  with lower id_num are
returned.


* **Parameters**

    **natoms** (*integer*) – The number of atoms to return.



* **Returns**

    **bonded_atoms** – The list with the index of natoms atoms bonded self.



* **Return type**

    list of int


<!-- !! processed by numpydoc !! -->

#### connect(self, atom: 'AtomTop')
Connects self with other atom setting the bond.


* **Parameters**

    **atom** (*AtomTop*) – Atom to connect.


<!-- !! processed by numpydoc !! -->

#### copy(self)
Returns a copy of the current atom.


* **Returns**

    **atom_top** – The copy of the atom.



* **Return type**

    AtomTop


<!-- !! processed by numpydoc !! -->

#### property residname()
An identifier of the residue (resid+name)


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



### class gaddlemaps.components.Atom(atom_top, atom_gro)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

An atom class that wraps the AtomTop and AtomGro classes.

You can access to the methods and attributes that both AtomGro and
AtomItp have. To create the atom object, both input atoms should have the
same resname and name attributes. On the other hand, only the attributes
from the AtomGro can be changed (e.g. positions, velocities, …) excluding
the resname and name.


* **Parameters**

    
    * **atom_top** (*AtomTop*) – The AtomTop object.


    * **atom_gro** (*AtomGro*) – The AtomGro object.



* **Raises**

    
    * [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If the atom_gro and atom_top do not correspond to the same atom.


    * [**TypeError**](https://docs.python.org/3/library/exceptions.html#TypeError) – If the inputs are not instances of the corresponding classes.



* **Attributes**

    `atom_gro`

        AtomGro : The input AtomGro object

    `atom_top`

        AtomTop : The input AtomTop object

    `gro_resid`

        Residue number for the gro part of the atom (atom_gro).

    `top_resid`

        Residue number for the part of the atom with the topology (atom_top).


### Methods

| `copy`(self)

                                  | Returns a copy of self.

                                                                                                                                                                     |
<!-- !! processed by numpydoc !! -->

#### property atom_gro()
The input AtomGro object


* **Type**

    AtomGro


<!-- !! processed by numpydoc !! -->

* **Return type**

    `AtomGro`



#### property atom_top()
The input AtomTop object


* **Type**

    AtomTop


<!-- !! processed by numpydoc !! -->

* **Return type**

    `AtomTop`



#### copy(self)
Returns a copy of self.

Only the gro atom is copied. The atom_top remains the same.


* **Returns**

    **new_atom** – The copied atom.



* **Return type**

    Atom


<!-- !! processed by numpydoc !! -->

#### property gro_resid()
Residue number for the gro part of the atom (atom_gro).

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`int`](https://docs.python.org/3/library/functions.html#int)



#### property top_resid()
Residue number for the part of the atom with the topology (atom_top).

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`int`](https://docs.python.org/3/library/functions.html#int)



### class gaddlemaps.components.Residue(atoms)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

A class with the information of a residue in a .gro file.

This class creates objects that are enumerations of AtomGro instances. It
has methods to manipulate atoms positions maintaining the shape of the
residue. This class have to be initialized with non empty list of atoms.

You can add two residues if the atoms from both residues have the same
residue name and number. This can also be done with just one atom and the
same check will be done.


* **Parameters**

    **atoms** (*list of AtomGro*) – A list with the atoms of the residue.


:raises ValueError : If the the input list of atoms is empty or they have different: residue number or name.


* **Attributes**

    `atoms`

        list of AtomGro: List with a copy of the atoms of the residue.

    `atoms_ids`

        list of int: A list with the ids of the atoms of the residue.

    `atoms_positions`

        numpy.ndarray((N, 3)) : An array with the atoms positions.

    `atoms_velocities`

        numpy.ndarray((N, 3)) or None : An array with the atoms velocities.

    `distance_to_zero`

        float : The distance between the geometric_center and (0, 0, 0)

    `geometric_center`

        numpy.ndarray(3): Coordinates of the geometric center of the residue.

    `resid`

        int: Residue number of the residue.

    `residname`

        string: An identifier of the residue (resid+name)

    `resname`

        string: Resname of the residue.

    `x`

        float: The x coordinate of the geometric center of the residue.

    `y`

        float: The y coordinate of the geometric center of the residue.

    `z`

        float: The z coordinate of the geometric center of the residue.


### Methods

| `copy`(self)

                                  | Returns a copy of the residue.

                                                                                                                                                              |
| `distance_to`(self, residue, numpy.ndarray], …)

 | Returns the distance between self and residue.

                                                                                                                                              |
| `move`(self, displacement)

                      | Moves the residue a given displacement vector.

                                                                                                                                              |
| `move_to`(self, new_position)

                   | Moves the residue geometric_center to new_position.

                                                                                                                                         |
| `remove_atom`(self, atom)

                       | Removes a given atom from the residue.

                                                                                                                                                      |
| `rotate`(self, rotation_matrix)

                 | Rotate the residue around its center of mass with a given rotation matrix.

                                                                                                                  |
| `update_from_molecule_top`(self, mtop)

          | Modifies the Residue atoms name to match the mtop.

                                                                                                                                          |
| `write_gro`(self, fout)

                         | Writes a .gro file with the residue conformation.

                                                                                                                                           |
<!-- !! processed by numpydoc !! -->

#### property atoms()
List with a copy of the atoms of the residue.


* **Type**

    list of AtomGro


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[`AtomGro`]



#### property atoms_ids()
A list with the ids of the atoms of the residue.


* **Type**

    list of int


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`int`](https://docs.python.org/3/library/functions.html#int)]



#### property atoms_positions()
An array with the atoms positions.


* **Type**

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)((N, 3))


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### property atoms_velocities()
An array with the atoms velocities.
If one of the atoms has no velocity this returns None.


* **Type**

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)((N, 3)) or [None](https://docs.python.org/3/library/constants.html#None)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)]



#### copy(self)
Returns a copy of the residue.


* **Returns**

    **residue** – The copy of the residue.



* **Return type**

    Residue


<!-- !! processed by numpydoc !! -->

#### distance_to(self, residue: Union[ForwardRef('Residue'), numpy.ndarray], box_vects: numpy.ndarray = None, inv: bool = False)
Returns the distance between self and residue.

residue can be a Residue instance or a 3D vector.


* **Parameters**

    
    * **residue** (*Residue** or *[*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)) – The residue or a point to compute the distance.


    * **box_vects** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)) – The box vectors to apply periodic boundary conditions.


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



#### property geometric_center()
Coordinates of the geometric center of the residue.


* **Type**

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)(3)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### move(self, displacement: numpy.ndarray)
Moves the residue a given displacement vector.


* **Parameters**

    **displacement** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**3**)*) – An array with the displacement vector.


<!-- !! processed by numpydoc !! -->

#### move_to(self, new_position: numpy.ndarray)
Moves the residue geometric_center to new_position.


* **Parameters**

    **new_position** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**3**)*) – An array with the new position coordinates.


<!-- !! processed by numpydoc !! -->

#### remove_atom(self, atom: 'AtomGro')
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



#### property resname()
Resname of the residue.


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### rotate(self, rotation_matrix: numpy.ndarray)
Rotate the residue around its center of mass with a given rotation
matrix.


* **Parameters**

    **rotation_matrix** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)) – 3x3 array with the rotation matrix.
    ADVICE: create the matrix with the help of “rotation_matrix”
    function placed in the root of the package.


<!-- !! processed by numpydoc !! -->

#### update_from_molecule_top(self, mtop: 'MoleculeTop')
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

#### write_gro(self, fout: str)
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



### class gaddlemaps.components.MoleculeTop(ftop, file_format=None)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Loads molecules from a topology file.

This class behaves like a list of atoms which has bonds defined. The
appropriate parser will be used based on the input file extension. The
available parsers are summarized in the class attribute “PARSERS”. In this
attribute, the keys are the files extensions and the values the
corresponding functions that extracts the information from the files with
that extensions. These functions should return:

> 
> * The name of the molecule


> * A list with tuples with the atoms and residues names in order of

> appearance in the file.
> - A list with tuples with atoms index (referred to the atoms_info
> indexes) that are bonded.


* **Parameters**

    
    * **ftop** (*string*) – The path to the file with the molecule name and bonds information.


    * **file_format** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)*, **Optional*) – The file extension of ftop. If it is None this will be taken from
    ftop.



* **Raises**

    
    * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If the file format is not supported.


    * [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If the input file misses information.



* **Attributes**

    `resids`

        list of str: Residue names of the atoms without consecutive

    `resname_len_list`

        list of tuple(str, int) :

    `resnames`

        list of str: Residue names of the atoms without consecutive


### Methods

| `copy`(self)

                                    | Returns a copy of the molecule_top.

                                                                                                                                                         |
| `index`(self, atom)

                             | Returns the index of the atom in the molecule.

                                                                                                                                              |
<!-- !! processed by numpydoc !! -->

#### copy(self)
Returns a copy of the molecule_top.

The atoms forming the copy are not the same objects as the original
molecule so you do not have to worry about linked objects.


* **Returns**

    **molecule_top** – The copy of the molecule.



* **Return type**

    MoleculeTop


<!-- !! processed by numpydoc !! -->

#### index(self, atom: 'AtomTop')
Returns the index of the atom in the molecule.


* **Parameters**

    **atom** (*AtomTop*) – The atom to find the index.



* **Returns**

    **index** – The index of the atom in the molecule.



* **Return type**

    [int](https://docs.python.org/3/library/functions.html#int)


<!-- !! processed by numpydoc !! -->

#### property resids()
Residue names of the atoms without consecutive

    repetitions.

To set this property a list with the same length of residues must be
passed.


* **Type**

    list of str


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`int`](https://docs.python.org/3/library/functions.html#int)]



#### property resname_len_list()
list of tuple(str, int):
[(resname_1, number_of_atoms_with_resname_1),
(resname_2, number_of_atoms_with_resname_2), …]

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int)]]



#### property resnames()
Residue names of the atoms without consecutive

    repetitions.

To set this property a list with the same length of residues must be
passed.


* **Type**

    list of str


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`str`](https://docs.python.org/3/library/stdtypes.html#str)]



### class gaddlemaps.components.Molecule(molecule_top, residues)
Bases: `gaddlemaps.components._residue.Residue`

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

| `copy`(self, new_residues)

                      | Returns a copy of the molecule.

                                                                                                                                                             |
| `deep_copy`(self, new_residues)

                 | Returns a deep copy of the molecule.

                                                                                                                                                        |
| `distance_to`(self, residue, numpy.ndarray], …)

 | Returns the distance between self and residue.

                                                                                                                                              |
| `from_files`(fgro, ftop)

                        | Loads the molecule from gro and a compatible topology file.

                                                                                                                                 |
| `index`(self, atom)

                             | Returns the index of the atom in the molecule.

                                                                                                                                              |
| `move`(self, displacement)

                      | Moves the residue a given displacement vector.

                                                                                                                                              |
| `move_to`(self, new_position)

                   | Moves the residue geometric_center to new_position.

                                                                                                                                         |
| `remove_atom`(self, atom)

                       | Removes a given atom from the residue.

                                                                                                                                                      |
| `rotate`(self, rotation_matrix)

                 | Rotate the residue around its center of mass with a given rotation matrix.

                                                                                                                  |
| `update_from_molecule_top`(self, mtop)

          | Modifies the Residue atoms name to match the mtop.

                                                                                                                                          |
| `write_gro`(self, fout)

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

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)((N, 3))


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### property atoms_velocities()
An array with the atoms velocities.
If one of the atoms has no velocity this returns None.


* **Type**

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)((N, 3)) or [None](https://docs.python.org/3/library/constants.html#None)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)]



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



#### copy(self, new_residues: List[gaddlemaps.components._residue.Residue] = None)
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

#### deep_copy(self, new_residues: List[gaddlemaps.components._residue.Residue] = None)
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

#### distance_to(self, residue: Union[ForwardRef('Residue'), numpy.ndarray], box_vects: numpy.ndarray = None, inv: bool = False)
Returns the distance between self and residue.

residue can be a Residue instance or a 3D vector.


* **Parameters**

    
    * **residue** (*Residue** or *[*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)) – The residue or a point to compute the distance.


    * **box_vects** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)) – The box vectors to apply periodic boundary conditions.


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



#### classmethod from_files(fgro: str, ftop: str)
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

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)(3)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### index(self, atom: 'Atom')
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



#### move(self, displacement: numpy.ndarray)
Moves the residue a given displacement vector.


* **Parameters**

    **displacement** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**3**)*) – An array with the displacement vector.


<!-- !! processed by numpydoc !! -->

#### move_to(self, new_position: numpy.ndarray)
Moves the residue geometric_center to new_position.


* **Parameters**

    **new_position** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**3**)*) – An array with the new position coordinates.


<!-- !! processed by numpydoc !! -->

#### remove_atom(self, atom: 'AtomGro')
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



#### rotate(self, rotation_matrix: numpy.ndarray)
Rotate the residue around its center of mass with a given rotation
matrix.


* **Parameters**

    **rotation_matrix** ([*numpy.ndarray*](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)) – 3x3 array with the rotation matrix.
    ADVICE: create the matrix with the help of “rotation_matrix”
    function placed in the root of the package.


<!-- !! processed by numpydoc !! -->

#### update_from_molecule_top(self, mtop: 'MoleculeTop')
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

#### write_gro(self, fout: str)
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



### class gaddlemaps.components.SystemGro(fgro)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Class to work with the information in gro files.

Basically this class acts like a list of Residue objects. Only one Residue
instance of each type is loaded to afford memory for large systems. The
positions of the rest of the residues are storage and they are generated
when requested.


* **Parameters**

    **fgro** (*string*) – Gromacs file with the system information.



* **Attributes**

    `box_matrix`

        numpy.ndarray (3,3) : the 3 lattice vectors of the box.

    `comment_line`

        str: The comment line in the gro file.

    `composition`

        Counter of str: int : For each resname (key), how many molecules there are (value).

    `molecules_info_ordered_all`

        generator of int: Returns the index of the molecules in the

    `molecules_resname_len_index`

        dict of tuple (string, int) to int: The keys of the dictionary are

    `n_atoms`

        int : The number of atoms in the system.


<!-- !! processed by numpydoc !! -->

#### property box_matrix()
the 3 lattice vectors of the box.


* **Type**

    [numpy.ndarray](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray) (3,3)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### property comment_line()
The comment line in the gro file.


* **Type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### property composition()
For each resname (key), how many molecules
there are (value).


* **Type**

    Counter of str



* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Mapping`](https://docs.python.org/3/library/typing.html#typing.Mapping)[[`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int)]



#### property molecules_info_ordered_all()
Returns the index of the molecules in the
different_molecules attribute in order of appearance.


* **Type**

    generator of int


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Generator`](https://docs.python.org/3/library/typing.html#typing.Generator)[[`int`](https://docs.python.org/3/library/functions.html#int), `None`, `None`]



#### property molecules_resname_len_index()
The keys of the dictionary are
tuples with the residue name and the number of atoms of the
different residues in the system and the values are its index in
the list of different molecules.


* **Type**

    dict of tuple (string, [int](https://docs.python.org/3/library/functions.html#int)) to int


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Dict`](https://docs.python.org/3/library/typing.html#typing.Dict)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int)], [`int`](https://docs.python.org/3/library/functions.html#int)]



#### property n_atoms()
The number of atoms in the system.


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`int`](https://docs.python.org/3/library/functions.html#int)



### class gaddlemaps.components.System(fgro, \*ftops)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

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

| `add_ftop`(self, ftop)

                          | Adds and identifies the molecule from the ftop to the system.

                                                                                                                               |
| `add_molecule_top`(self, mol_top)

               | Adds a molecule to the system and find it in the gro file.

                                                                                                                                  |
<!-- !! processed by numpydoc !! -->

#### add_ftop(self, ftop: str)
Adds and identifies the molecule from the ftop to the system.

<!-- !! processed by numpydoc !! -->

#### add_molecule_top(self, mol_top: gaddlemaps.components._components_top.MoleculeTop)
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



### gaddlemaps.components.are_connected(atoms: Union[List[gaddlemaps.components._components.Atom], List[gaddlemaps.components._components_top.AtomTop]])
Check if the input atoms are connected.


* **Parameters**

    **atoms** (*List of AtomTop** or **Atom*) – The list of atoms to compute if they are connected.



* **Returns**

    **connected** – True if the atoms are connected, else False.



* **Return type**

    [bool](https://docs.python.org/3/library/functions.html#bool)


<!-- !! processed by numpydoc !! -->
