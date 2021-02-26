# gaddlemaps package

## Subpackages


* gaddlemaps.components package


    * Module contents


        * Molecular Simulation Components


* gaddlemaps.parsers package


    * Module contents


        * Simulation Files Parsers


## Module contents

### gaddlemaps — Change molecules in simulations with python

This python package provides an implementation of the GADDLE-Maps (General
Algorithm for Discrete Object Deformations Based on Local Exchange Maps) that
allows to change molecules in a molecular dynamics simulations. For example,
this tool is very handy to back-map a coarse grained simulation to atomistic
forcefields and also in the other way around.

TODO: Complete this description

<!-- !! processed by numpydoc !! -->

### class gaddlemaps.Alignment(start=None, end=None)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Class to manage the molecules alignment.

This class serve as interface during the molecular alignment process. It
will take two molecules in different representations (e.g. coarse grained
and atomistic), one of them will be taken as initial representation and the
other as the final one. These molecules can be passed as arguments in the
object initialization or assigned after its creation.

NOTE: Once both molecules are assigned you can only set them again with the
same molecule type (although it may have different position).

Once the molecules are fixed you can call the “aling_molecules” method to
find the optimum overlap between both representation. Then, you should
call the “init_exchange_map” method to generate the map to change between
resolution for other start molecules configurations. You can also call
the “write_comparative_gro” method to write a .gro file with the two
aligned molecules with different residue names to visualize the found
overlap in other visualization programs.


* **Parameters**

    
    * **start** (*Molecule*) – The molecule to align in the initial resolution.


    * **end** (*Molecule*) – The molecule to align in the final resolution.



#### exchange_map()
The map to change from initial resolution to the end one. It has to be
initialized with the init_exchange_map method.


* **Type**

    ExchangeMap



#### SIGMA_SCALE()
A number that modulates the molecule displacements.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)



#### STEPS_FACTOR()
The factor to apply in the calculation of the number of steps in the
Monte-Carlo alignment.


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)



* **Attributes**

    `end`

        Molecule : The molecule in the final resolution.

    `start`

        Molecule : The molecule in the initial resolution.


### Methods

| `align_molecules`([restrictions, …])

 | Starts the alignment engine to find the optimal overlap between molecs.

 |
| `init_exchange_map`([scale_factor])

  | Initializes the exchange map with the current molecules configuration.

  |
| `interactive_restrictions`()

         | Creates the widget to visually generate the restrictions for alignment.

 |
| `write_comparative_gro`([fname])

     | Writes a .gro file with start and end molecules to check the overlap.

   |
<!-- !! processed by numpydoc !! -->

#### SIGMA_SCALE( = 0.5)

#### STEPS_FACTOR( = 5000)

#### align_molecules(restrictions=None, deformation_types=None, ignore_hydrogens=True, auto_guess_protein_restrictions=True)
Starts the alignment engine to find the optimal overlap between molecs.

NOTE: If None is input as restrictions and the molecules to align
have multiple residues and also the auto_guess_protein_restrictions
parameter is True, they will be guessed by residue matching (see
guess_protein_restrains function). However take into account that the
time taken by this step scales with the square of the number of atoms if
no restrictions are given and linearly if they are guessed for all the
atoms by residue matching.

If you know that the molecules in both resolution are in the same
configuration you may want to perform the alignment just rotating and
translating the molecules as a whole and avoid molecular deformation by
setting the “deformation_types” parameter to (0, 1). This will translate
in a better performance.

On the other hand, in the most of the cases, the hydrogen atoms are not
very relevant atoms for the alignment so they are by default (see
“ignore_hydrogens” parameter) not taken into account in the distance
minimization process when they belong to the molecule with larger number
of atoms (the one that will remain still in the process).


* **Parameters**

    
    * **restrictions** (*list of tuple of int**, **optional*) – A list of tuples with pairs of atom numbers corresponding to
    start and end atoms indexes in the molecules. The align will be
    performed privileging configurations where those atoms are close.
    By default is set to [].

    Example:
    >>> restrictions = [(1, 3), (4, 5)]

    **IMPORTANT INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
    (STARTS IN 0).**



    * **deformation_types** (*tuple of int**, **optional*) – Specifies the type of the minimization. Possible options:

        0 : Translation
        1 : Rotation
        2 : Individual atom move

    If it is None, all the possibilities are chosen.



    * **ignore_hydrogens** ([*bool*](https://docs.python.org/3/library/functions.html#bool)*, **optional*) – If True, hydrogen atoms will not be included in the minimization
    of the distances. This will be only applied to the molecule which
    is not moved in the alignment engine.


    * **auto_guess_protein_restrictions** ([*bool*](https://docs.python.org/3/library/functions.html#bool)*, **optional*) – If True automatic restrictions will try to be guessed if the
    molecules to be aligned have multiple residues.



* **Raises**

    [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If the automatic restrictions guess failed. In that case, consider
    to set “auto_guess_protein_restrictions” parameter to False.


<!-- !! processed by numpydoc !! -->

#### property end()
The molecule in the final resolution.

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[`Molecule`]



#### init_exchange_map(scale_factor=0.5)
Initializes the exchange map with the current molecules configuration.


* **Parameters**

    **scale_factor** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The compression factor to apply to mapped molecules (see the
    ExchangeMap class documentation for more information).


<!-- !! processed by numpydoc !! -->

#### interactive_restrictions()
Creates the widget to visually generate the restrictions for alignment.


* **Return type**

    [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`ForwardRef`](https://docs.python.org/3/library/typing.html#typing.ForwardRef), [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`int`](https://docs.python.org/3/library/functions.html#int)]]]]



* **Returns**

    
    * **restriction_widget** (*ipywidgets.Widget*) – The widget that contains the constraint generator.


    * **restrictions** (*List[Tuple[int, int]]*) – The restrictions that will be generated by the widget. This will
    be initially empty and it will be filled as the widget is used.



<!-- !! processed by numpydoc !! -->

#### property start()
The molecule in the initial resolution.

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[`Molecule`]



#### write_comparative_gro(fname=None)
Writes a .gro file with start and end molecules to check the overlap.

The molecules have START and END residue names respectively to ease the
representation.


* **Parameters**

    **fname** (*string**, **optional*) – The output .gro file name. If it is None, fname is set as:
    {name}_compare.gro, where name is the name of the molecule.


<!-- !! processed by numpydoc !! -->

### class gaddlemaps.Chi2Calculator(mol1, mol2, restrictions=None)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Functor that calculates the chi2 between 2 positions arrays.

This class is initialized with the restrictions that must be taken into
account. This avoid repeating array accessing. Once the object is
initialized it can be called with new atomic positions  to calculate the new
chi2 value associated with the distance between set of points. This distance
is calculated quite differently depending on the given restrictions. In case
of no restrictions see “chi2_molecules” method. If some restrictions are
given the distances between atoms is calculated between of pairs of atoms
specified in the restrictions instead of between closest atoms.


* **Parameters**

    
    * **mol1** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**(**N**, **3**)**)*) – Array with the atomic positions of the molecule that will be still.
    These positions will remain constant in  future object call.


    * **mol2** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**(**N**, **3**)**)*) – Array with the atomic positions of the molecule that will be changing
    during the optimum overlap finding process.


    * **restriction** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**(**N**, **2**)**) or **array convertible**, **optional*) – A list of tuples with pairs of atom numbers corresponding to mol1
    and mol2 atoms molecules. The align will be performed privileging
    configurations where those atoms are close. By default is set to [].

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE**



### Methods

| `__call__`(mol2)

                     | Call self as a function.

                                                |
| `chi2_molecules`(mol2)

               | Computes the chi2 by calculating the distance between molecules.

        |
<!-- !! processed by numpydoc !! -->

#### chi2_molecules(mol2)
Computes the chi2 by calculating the distance between molecules.

This function computes a chi2 based on the distance between nearest
atoms of two molecules. Basically sums the distance between the atoms
of mol2 and the nearest atom of mol1.


* **Parameters**

    **mol2** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**(**N**,**3**)**)*) – Position of the second molecule atoms.



* **Returns**

    **chi2** – Computed chi2.



* **Return type**

    [float](https://docs.python.org/3/library/functions.html#float)


<!-- !! processed by numpydoc !! -->

### class gaddlemaps.ExchangeMap(refmolecule, targetmolecule, scale_factor=0.5)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Functor to extrapolate atomic resolution to other configurations.

When this class is initialized, a functor is created. It has to be
initialized with the molecules in the initial and final resolution
overlapped. Then you can call this method with new molecules in the
initial resolution to obtain its representation in the final one. The new
coordinates of the extrapolated molecules are scaled by the “scale_factor”.
This factor should be smaller than one if you want to extrapolate a complete
system for future simulations. This avoids molecular overlapping and
prevents the simulations to crash.


* **Parameters**

    
    * **refmolecule** (*Molecule*) – Molecule in the initial resolution.


    * **targetmolecule** (*Molecule*) – Molecule in the final resolution.


    * **scale_factor** ([*float*](https://docs.python.org/3/library/functions.html#float)*, **Optional*) – The factor that modulates the scale of the atoms positions in the new
    resolution respect it closest atom in the initial resolution.


### Examples

```python
>>> BmimCG_ref = Molecule(fgro, fitp)
>>> BmimAA_target = Molecule(groAA, itpAA)
```

```python
>>> transformation = ExchangeMap(BmimCG_ref, BmimAA_target)
```

```python
>>> BmimCG_new = Molecule(fgro2, fitp)
>>> BmimAA_new = transformation(BmimCG_new)
```


* **Attributes**

    `equivalences`

        dict of int to list of int : {r1_atom_index: [closest_r2_atoms_indexs]}


### Methods

| `__call__`(refmolecule)

              | This function takes as argument a molecule like the refmolecule, but in other position and returns its targetmolecule equivalent as a new molecule with the same res_number than the input.

 |
<!-- !! processed by numpydoc !! -->

#### property equivalences()
[closest_r2_atoms_indexs]}


* **Type**

    dict of int to list of int



* **Type**

    {r1_atom_index


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Dict`](https://docs.python.org/3/library/typing.html#typing.Dict)[[`int`](https://docs.python.org/3/library/functions.html#int), [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`int`](https://docs.python.org/3/library/functions.html#int)]]



### class gaddlemaps.Manager(system)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

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

### gaddlemaps.accept_metropolis(energy_0, energy_1, acceptance=0.01)
Evaluates if the new configuration with energy_1 is accepted or not.

This function returns True if the configuration with energy_1 is accepted
and False otherwise. If energy_1 <= energy_0 the configuration is accepted.
Else, it will be also accepted if a random number between 0 and 1 is
smaller than acceptance\*energy_0/energy_1.


* **Parameters**

    
    * **energy_0** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The value of the energy before the change.


    * **energy_1** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The value of the energy after the change.


    * **acceptance** ([*float*](https://docs.python.org/3/library/functions.html#float)* (**Optional**)*) – acceptance factor that regulates how many unfavorable cases are accepted



* **Returns**

    **change** – True if the change is accepted and False if not.



* **Return type**

    Bool


<!-- !! processed by numpydoc !! -->

### gaddlemaps.calcule_base(pos)
Calculates a orthonormal base from a list of 3 atoms.

Given a list of three vectors with the position of three atoms, this
function returns a vector basis and the application point. The first
vector goes from the application point to the last atom. The second one
is normal to the plane which contains the three atoms and perpendicular
to the first vector. The last is perpendicular to the others. All are
unitary forming an orthonormal basis. In case of collinear points, the
second vector is set to ([vec1[1], -vec1[0], 0]).


* **Parameters**

    **pos** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – List with three vectors (numpy.ndarray) with the position of the
    atoms.



* **Return type**

    [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), …], [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)]



* **Returns**

    
    * **base** (*tuple of numpy.ndarray*) – A tuple with the three vectors of the base.


    * **app_point** (*float*) – The application point of the base. It corresponds to pos[0]



<!-- !! processed by numpydoc !! -->

### gaddlemaps.find_atom_random_displ(atoms_pos, bonds_info, atom_index, sigma_scale=0.5)
Finds a random displacement for the atom with a given index.

This displacement is chosen in a perpendicular direction according to the
number of bonded atoms.


* **Parameters**

    
    * **atoms_pos** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – An array with the positions of the molecule atoms in rows.


    * **bonds_info** (*dictionary*) – A dict with the information of the bonds. Example:

        bonds_info = {0:[(1, 5.4), (2, 6.4)], 1:[(0, 5.4), ], …}

    The keys refers to atom index and values are lists with tuples. Each
    tuple contains the bonded atom index and the bond length.



    * **atom_index** (*integer** (**Optional**)*) – The index of the atom to move (respecting the index of atoms_pos). If
    None is given a random one is taken.


    * **sigma_scale** ([*float*](https://docs.python.org/3/library/functions.html#float)) – A factor to scale the sigma of the distribution of the displacement
    module.



* **Returns**

    **displ** – The displacement vector to sum to the position of the interest atom.



* **Return type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)


<!-- !! processed by numpydoc !! -->

### gaddlemaps.guess_protein_restrains(mol1, mol2)
Guess restrains for molecules with multiple residues.

Checks if mol1 and mol2 have the same residue names and create restrains
pairing atoms with the same residues. This function only works with
molecules with multiple residues.


* **Parameters**

    
    * **mol1** (*Molecule** or **MoleculeGro*) – The first molecule to find the restrains.


    * **mol2** (*Molecule** or **MoleculeGro*) – The second molecule to find the restrains.



* **Returns**

    **restrains** – A list of tuples with pairs of atom numbers corresponding to start
    and end atoms molecules. The align will be performed privileging
    configurations where those atoms are close. By default is set to [].

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
    (STARTS IN 0).**




* **Return type**

    list of tuple of int



* **Raises**

    [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If the input molecules have different residue names.


<!-- !! processed by numpydoc !! -->

### gaddlemaps.guess_residue_restrains(res1, res2, offset1=0, offset2=0)
Guess restrains for molecules with one residue.


* **Parameters**

    
    * **res1** (*Residue*) – The first residue to find the restrains.


    * **res2** (*Residue*) – The second residue to find the restrains.


    * **offset1** ([*int*](https://docs.python.org/3/library/functions.html#int)) – An offset to add to the atom index of res1.


    * **offset2** ([*int*](https://docs.python.org/3/library/functions.html#int)) – An offset to add to the atom index of res2.



* **Returns**

    **restrains** – A list of tuples with pairs of atom numbers corresponding to start
    and end atoms residues. The align will be performed privileging
    configurations where those atoms are close. By default is set to [].

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
    (STARTS IN 0).**




* **Return type**

    list of tuple of int


<!-- !! processed by numpydoc !! -->

### gaddlemaps.interactive_restrictions(correspondence, style=None)
Creates the widget to generate the restrictions of all the species in the
alignment. It generates the final representation fo the widget.


* **Parameters**

    
    * **manager** (*Manager*) – The object that manages all the alignment process


    * **style** (*Optional**[*[*int*](https://docs.python.org/3/library/functions.html#int)*]*) – An integer that determine which style will be used to represent the
    widget for each specie.

    > 0: One tab per specie.
    > 1: Accordion, when one specie opens the other collapse
    > 2: Vertically aligned, one over the other

    The default value is 2. This is the only one fully operational, in the
    other ones it is necessary to manually refresh the widget in the
    notebook when changing between species.




* **Return type**

    [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[`Widget`, [`Dict`](https://docs.python.org/3/library/typing.html#typing.Dict)[[`str`](https://docs.python.org/3/library/stdtypes.html#str), [`Optional`](https://docs.python.org/3/library/typing.html#typing.Optional)[[`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`int`](https://docs.python.org/3/library/functions.html#int)]]]]]



* **Returns**

    
    * **restriction_widget** (*ipywidgets.Widget*) – The widget that contains the constraint generator for all the species


    * **restrictions** (*Dict[str, List[Tuple[int, int]]]*) – The dictionary with the restrictions that will be generated by the
    widget for each specie. This will be initially empty and it will be
    filled as the widget is used.



<!-- !! processed by numpydoc !! -->

### gaddlemaps.minimize_molecules(mol1_positions, mol2_positions, mol2_com, sigma_scale, n_steps, restriction, mol2_bonds_info, displacement_module, sim_type)
Minimizes the distance between two molecules.


* **Parameters**

    
    * **mol1_positions** (*np.array**(**(**N**, **3**)**)*) – The positions of the atoms of the static molecule.


    * **mol2_positions** (*np.array**(**(**N**, **3**)**)*) – The positions of the atoms of the mobile molecule.


    * **mol2_com** (*np.array**(**3**)*) – The center of mass (or the geometric center) of the mobile molecule.
    This will be used to apply de rotations.


    * **sigma_scale** ([*float*](https://docs.python.org/3/library/functions.html#float)) – A number that module the amplitude of the single atom displacements.


    * **n_steps** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The number of steps without changes in the chi2 in the Monte-Carlo
    alignment.


    * **restriction** (*list of tuple of int*) – A list of tuples with pairs of atom numbers corresponding to mol1
    and mol2 atoms molecules. The align will be performed privileging
    configurations where those atoms are close.

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE
    (STARTS IN 0).**



    * **mol2_bonds_info** (*defaultdict of int: list of tuple**(*[*int*](https://docs.python.org/3/library/functions.html#int)*, *[*float*](https://docs.python.org/3/library/functions.html#float)*)*) – A complex data structure that collect the information of the bond
    distances. The key of the property corresponds to the atom index in the
    molecule. The value is a list with tuples. For each tuple, the first
    value corresponds with the index of the bonded atom and the second is
    the length of the bond.


    * **same_com** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, translations are not allowed.


    * **anchura** ([*float*](https://docs.python.org/3/library/functions.html#float)) – A number that modules the amplitude in the translations displacements.


    * **sim_type** (*tuple of int*) – Specifies the type of the minimization. Possible options:

        0 : Translation
        1 : Rotation
        2 : Individual atom move

    If it is a list, all specified methods are combined.



<!-- !! processed by numpydoc !! -->

### gaddlemaps.move_mol_atom(atoms_pos, bonds_info, atom_index=None, displ=None, sigma_scale=0.5)
Moves an atom of a molecule respecting almost all bond distances.

By default, a random atom is picked from atom_pos and moved randomly in a
certain direction (see the published article for a better description of
this step).


* **Parameters**

    
    * **atoms_pos** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – An array with the positions of the molecule atoms in rows.


    * **bonds_info** (*dictionary*) – A dict with the information of the bonds. Example:

        bonds_info = {0:[(1, 5.4), (2, 6.4)], 1:[(0, 5.4), ], …}

    The keys refers to atom index and values are lists with tuples. Each
    tuple contains the bonded atom index and the bond length.



    * **atom_index** (*integer** (**Optional**)*) – The index of the atom to move (respecting the index of atoms_pos). If
    None is given a random one is taken.


    * **displ** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)* (**Optional**)*) – The displacement vector. If None is given a random displacement is
    calculated in the normal plane to the line jointing most nearest atoms.


    * **sigma_scale** ([*float*](https://docs.python.org/3/library/functions.html#float)) – A factor to scale the sigma of the distribution of the
    displacement module.



* **Returns**

    **modified_atoms_pos** – An array with the modified positions of the atoms in rows.



* **Return type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)


<!-- !! processed by numpydoc !! -->

### gaddlemaps.remove_hydrogens(molecule, restrictions)
Returns positions of atoms that are not hydrogens and fixes restrictions.


* **Parameters**

    
    * **molecule** (*Molecule*) – The molecule to remove hydrogens.


    * **restrictions** (*list of tuple of int*) – A list of tuples with pairs of atom numbers corresponding to start
    and end atoms molecules. The align will be performed privileging
    configurations where those atoms are close. By default is set to [].

    ### Example

    ```python
    >>> restrictions = [(1, 3), (4, 5)]
    ```

    **IMPORTANT: INDEX ARE REFERENCED TO THE ATOM INDEX IN THE MOLECULE**




* **Return type**

    [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), [`List`](https://docs.python.org/3/library/typing.html#typing.List)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`int`](https://docs.python.org/3/library/functions.html#int)]]]



* **Returns**

    
    * **new_positions** (*numpy.ndarray*) – The positions of the atoms that are not hydrogens.


    * **new_restrictions** (*list of tuple of int*) – The restrictions with the correct indexes.



<!-- !! processed by numpydoc !! -->

### gaddlemaps.rotation_matrix(axis, theta)
Returns the rotation matrix associated to an angle and an axis.


* **Parameters**

    
    * **axis** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A vector in 3D space that defines the axis of rotation


    * **theta** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The angle to rotate in radians in counter clockwise.


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)
