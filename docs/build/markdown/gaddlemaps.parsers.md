# gaddlemaps.parsers package

## Module contents

### Simulation Files Parsers

This submodule contains some useful functionalities to read and extract
information from files coming from molecular simulation programs. Up to now,
only files in the format used by GROMACS are contemplated (.gro and .itp
formats).

<!-- !! processed by numpydoc !! -->

### class gaddlemaps.parsers.CoordinatesParser(path, mode='r')
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)


* **Attributes**

    **EXTENSIONS**

    `box_matrix`

        Return a 3x3 matrix with the 3 lattice vectors

    `comment`

        A comment/name of the system

    `natoms`

        Number of atoms in the file


### Methods

| `close`()

 | Closes the file

 |
| `next`()

                                  | Returns next atomline formatted

                                                                                                                                                             |
| `seek_atom`(index)

                        | Displaces the position of the ‘cursor’ to an atom line

                                                                                                                                      |
| `writeline`(atomlist)

                     | Writes a line of atom information

                                                                                                                                                           |
| `writelines`(list_atomlist)

               | Writes several lines of atoms

                                                                                                                                                               |
<!-- !! processed by numpydoc !! -->

#### EXTENSIONS(: Optional[Tuple[[str](https://docs.python.org/3/library/stdtypes.html#str), …]] = None)

#### abstract property box_matrix()
Return a 3x3 matrix with the 3 lattice vectors

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### abstract close()
Closes the file

<!-- !! processed by numpydoc !! -->

#### abstract property comment()
A comment/name of the system

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### abstract property natoms()
Number of atoms in the file

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`int`](https://docs.python.org/3/library/functions.html#int)



#### abstract next()
Returns next atomline formatted

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Union`](https://docs.python.org/3/library/typing.html#typing.Union)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float)], [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float)]]



#### abstract seek_atom(index)
Displaces the position of the ‘cursor’ to an atom line

Displaces the position of the ‘cursor’ to the beginning of the
line of the ‘index’ atom, where the first atom index is 0. This function
warrants that calling the “next” function the information about atom
“index” will be returned.


* **Parameters**

    **index** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The index of the atom to found


<!-- !! processed by numpydoc !! -->

#### abstract writeline(atomlist)
Writes a line of atom information

If there was no content written in the file it creates the
header and the number of atoms. If the number of atoms was not
provided it will kept empty and the number will be written just
before closing the file. This possibility is only compatible
with a number of atoms smaller than 1000000.


* **Parameters**

    **atomlist** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or **string*) – A list with all the info of the atom, just like the one returned by
    readline. If it is a string, it will be written directly without
    parsing.


<!-- !! processed by numpydoc !! -->

#### writelines(list_atomlist)
Writes several lines of atoms


* **Parameters**

    **list_atomlist** (*list of list of str** or *[*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – A list of lists with all the info of the atom, just like
    the one returned by readline


<!-- !! processed by numpydoc !! -->

### class gaddlemaps.parsers.GroFile(path, mode='r')
Bases: `gaddlemaps.parsers.CoordinatesParser`

Implements a file object for opening and writing .gro files

Grofile verifies the gromacs format, and autodetects the correct
way to read the file. It modifies the methods read and write of a
common file with methods that take as input lists with all the info
of an atom. The readlines and writelines are modified to take
(return) lists with the info of an atom in lists (as returned by
the read function).

This class can also be initiated with an already opened file.

**NOTE**: The “+” modifier for the open mode is not allowed right now


* **Parameters**

    
    * **path** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The path to the file to be opened.


    * **mode** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – “r” to open the file in read mode and “w” to open it in write mode.



* **Raises**

    [**IOError**](https://docs.python.org/3/library/exceptions.html#IOError) – If the gro file has not the correct format



* **Attributes**

    `box_matrix`

        numpy.ndarray(3,3): the 3 lattice vectors of the box.

    `comment`

        string: The comment of the .gro file.

    `name`

        str: The filename of the open file

    `natoms`

        integer: The number of atoms in the system.

    `position_format`

        tuple of int : the number of total figures and decimal places for the


### Methods

| `close`()

                                 | Closes the file

                                                                                                                                                                             |
| `determine_format`(atomline)

              | Returns the format of the postion coordinates and whether or not the .gro file has velocities

                                                                                               |
| `next`()

                                  | Returns next atomline formatted

                                                                                                                                                             |
| `parse_atomline`(atomline[, format_dict])

 | Parses an atom line and returns its content in a list.

                                                                                                                                      |
| `parse_atomlist`(atomlist[, format_dict])

 | Convert a list of atom info to string with the appropriate format

                                                                                                                           |
| `readline`([parsed])

                      | Returns the next line of the gro file.

                                                                                                                                                      |
| `readlines`()

                             | Returns all the lines of a grofile parsed in lists.

                                                                                                                                         |
| `seek_atom`(index)

                        | Displaces the position of the ‘cursor’ to an atom line

                                                                                                                                      |
| `validate_string`(string)

                 | Validates a string to be valid as resname or as name.

                                                                                                                                       |
| `writeline`(atomlist)

                     | Writes a line of atom information

                                                                                                                                                           |
| `writelines`(list_atomlist)

               | Writes several lines of atoms

                                                                                                                                                               |
<!-- !! processed by numpydoc !! -->

#### COORD_START( = 20)

#### DEFAULT_COMMENT( = "Gro file genereted with 'Gromacs Tools' python module.")

#### DEFAULT_POSTION_FORMAT( = (8, 3))

#### EXTENSIONS(: Optional[Tuple[[str](https://docs.python.org/3/library/stdtypes.html#str), …]] = ('gro', 'GRO'))

#### NUMBER_FIGURES( = 9)

#### property box_matrix()
the 3 lattice vectors of the box.


* **Type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)(3,3)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### close()
Closes the file

If the file is in write mode it will write al the remaining
information.

<!-- !! processed by numpydoc !! -->

#### property comment()
The comment of the .gro file.


* **Type**

    string


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### classmethod determine_format(atomline)
Returns the format of the postion coordinates and whether or
not the .gro file has velocities


* **Parameters**

    **atomline** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – An atom line from the .gro file



* **Returns**

    **form** – Dictionary with 2 keys, “position” and “velocities”. The
    position key leads to a tuple with 2 ints. THe first in is
    the number of characters in the float format (C) and the
    second is the number of decimal places in the float (D),
    i.e. corresponds to the ‘%C.Df’ format.



* **Return type**

    [dict](https://docs.python.org/3/library/stdtypes.html#dict)


<!-- !! processed by numpydoc !! -->

#### property name()
The filename of the open file


* **Type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`str`](https://docs.python.org/3/library/stdtypes.html#str)



#### property natoms()
The number of atoms in the system.


* **Type**

    integer


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`int`](https://docs.python.org/3/library/functions.html#int)



#### next()
Returns next atomline formatted

<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Union`](https://docs.python.org/3/library/typing.html#typing.Union)[[`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float)], [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`str`](https://docs.python.org/3/library/stdtypes.html#str), [`int`](https://docs.python.org/3/library/functions.html#int), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float), [`float`](https://docs.python.org/3/library/functions.html#float)]]



#### classmethod parse_atomline(atomline, format_dict=None)
Parses an atom line and returns its content in a list.

Parses an atom line and returns its content in a list. In its default
behaviour the input line is analyzed and the number of decimals are
guessed. It is also possible to insert the format in a format
dictionary. This is a dictionary with 2 keys, “position” and
“velocities”. The position key leads to a tuple with 2 ints. The
first in is the number of characters in the float format (C) and the
second is the number of decimal places in the float (D), i.e.
corresponds to the ‘%C.Df’ format.


* **Parameters**

    
    * **atomline** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – An string with a line from a grofile


    * **format_dict** (*dict_format** (**optional**)*) – None or dictionary with 2 keys, “position” and “velocities”. The
    position key leads to a tuple with 2 ints. THe first in is
    the number of characters in the float format (C) and the
    second is the number of decimal places in the float (D),
    i.e. corresponds to the ‘%C.Df’ format. If None the format
    will be automatically guessed.
    Default None.



* **Returns**

    **info** –

    A list with the following information:


        * mol_index(integer)


        * resname (string)


        * name (string)


        * global_index (integer)


        * x, y, z (floats)


        * vx, vy ,vz (floats) [optional]




* **Return type**

    [list](https://docs.python.org/3/library/stdtypes.html#list)


<!-- !! processed by numpydoc !! -->

#### classmethod parse_atomlist(atomlist, format_dict=None)
Convert a list of atom info to string with the appropriate format

Parses a list with the atom information and returns its content in a
string with the correct format to write in a .gro file. In its default
behaviour the input line is analyzed and the number of decimals are
guessed. It is also possible to insert the format in a format
dictionary. This is a dictionary with 2 keys, “position” and
“velocities”. The position key leads to a tuple with 2 ints. The
first is the number of characters in the float format (C) and the
second is the number of decimal places in the float (D), i.e.
corresponds to the ‘%C.Df’ format.


* **Parameters**

    
    * **atomlist** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – A list with the following information:


            * mol_index(integer)


            * resname (string)


            * name (string)


            * global_index (integer)


            * x, y, z (floats)


            * vx, vy ,vz (floats) [optional]]



    * **format_dict** ([*dict*](https://docs.python.org/3/library/stdtypes.html#dict)* (**optional**)*) – None or dictionary with 2 keys, “position” and “velocities”. The
    position key leads to a tuple with 2 ints. THe first in is
    the number of characters in the float format (C) and the
    second is the number of decimal places in the float (D),
    i.e. corresponds to the ‘%C.Df’ format. If None the format
    will be automatically guessed.
    Default None.



* **Returns**

    **atomline** – An string with a line from a .gro file



* **Return type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)


<!-- !! processed by numpydoc !! -->

#### property position_format()
the number of total figures and decimal places for the
positions.


* **Type**

    tuple of int


<!-- !! processed by numpydoc !! -->

* **Return type**

    [`Tuple`](https://docs.python.org/3/library/typing.html#typing.Tuple)[[`int`](https://docs.python.org/3/library/functions.html#int), [`int`](https://docs.python.org/3/library/functions.html#int)]



#### readline(parsed=True)
Returns the next line of the gro file.

If parsed=True (default) the data on the atom line is parsed
and returned in a list. If false the line is returned as an
string (just like if it were opened as a common file.)


* **Parameters**

    **parsed** ([*bool*](https://docs.python.org/3/library/functions.html#bool)* (**optional**)*) – Wether the data of the line should be parsed or not.
    Default: True



* **Returns**

    **info** – Returns the info of the line, i alist if it was parsed or
    in a string if not



* **Return type**

    [str](https://docs.python.org/3/library/stdtypes.html#str) or [list](https://docs.python.org/3/library/stdtypes.html#list)


<!-- !! processed by numpydoc !! -->

#### readlines()
Returns all the lines of a grofile parsed in lists.


* **Returns**

    **info** – Returns a list with the info of all the atom lines in lists



* **Return type**

    [list](https://docs.python.org/3/library/stdtypes.html#list)([list](https://docs.python.org/3/library/stdtypes.html#list))


<!-- !! processed by numpydoc !! -->

#### seek_atom(index)
Displaces the position of the ‘cursor’ to an atom line

Displaces the position of the ‘cursor’ to the beginning of the
line of the ‘index’ atom, where the first atom index is 0. If
the index is equal to the number of atoms the beginning of the
box lattice line is found.


* **Parameters**

    **index** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The index of the atom to found


<!-- !! processed by numpydoc !! -->

#### static validate_string(string)
Validates a string to be valid as resname or as name.

If the input name has more than 5 characters, new name is returned
cutting the input one.


* **Parameters**

    **string** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The string to validate.



* **Returns**

    **new_string** – The properly formatted string.



* **Return type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)


<!-- !! processed by numpydoc !! -->

#### writeline(atomlist)
Writes a line of atom information

If there was no content written in the file it creates the
header and the number of atoms. If the number of atoms was not
provided it will kept empty and the number will be written just
before closing the file. This possibility is only compatible
with a number of atoms smaller than 1000000.


* **Parameters**

    **atomlist** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or **string*) – A list with all the info of the atom, just like the one returned by
    readline. If it is a string, it will be written directly without
    parsing.


<!-- !! processed by numpydoc !! -->

#### writelines(list_atomlist)
Writes several lines of atoms

If there was no content written in the file it creates the
header and the number of atoms. If the number of atoms was not
provided it will kept empty and the number will be written just
before closing the file. This possibility is only compatible
with a number of atoms smaller than 1000000.


* **Parameters**

    **list_atomlist** (*list of list of str** or *[*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – A list of lists with all the info of the atom, just like
    the one returned by readline


<!-- !! processed by numpydoc !! -->

### class gaddlemaps.parsers.ParserManager()
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

### Methods

| **register**

                                |                                                                                                                                                                                             |
<!-- !! processed by numpydoc !! -->

#### parsers(: Dict[[str](https://docs.python.org/3/library/stdtypes.html#str), Type[gaddlemaps.parsers.CoordinatesParser]] = {'GRO': <class 'gaddlemaps.parsers.GroFile'>, 'gro': <class 'gaddlemaps.parsers.GroFile'>})

#### classmethod register(parser)
<!-- !! processed by numpydoc !! -->

### class gaddlemaps.parsers.ParserRegistered(name, bases, attrs)
Bases: [`abc.ABCMeta`](https://docs.python.org/3/library/abc.html#abc.ABCMeta)

### Methods

| `__call__`(\*args, \*\*kwargs)

               | Call self as a function.

                                                                                                                                                                    |
| `mro`(/)

                                  | Return a type’s method resolution order.

                                                                                                                                                    |
| `register`(subclass)

                      | Register a virtual subclass of an ABC.

                                                                                                                                                      |
<!-- !! processed by numpydoc !! -->

#### mro(/)
Return a type’s method resolution order.

<!-- !! processed by numpydoc !! -->

#### register(subclass)
Register a virtual subclass of an ABC.

Returns the subclass, to allow usage as a class decorator.

<!-- !! processed by numpydoc !! -->

### gaddlemaps.parsers.dump_lattice_gro(vectors)
Extracts final line of a .gro file from the lattice vectors for pbc.


* **Parameters**

    **vectors** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)*(**(**3**,**3**)**)*) – An array with a lattice vector in each row.



* **Returns**

    **line** – The final line of a .gro file.



* **Return type**

    string


<!-- !! processed by numpydoc !! -->

### gaddlemaps.parsers.extract_lattice_gro(line)
Extracts the lattice vectors for pbc from the final line of a .gro file.


* **Parameters**

    **line** (*string*) – The final line of a .gro file.



* **Returns**

    **vectors** – An array with a lattice vector in each row.



* **Return type**

    [numpy.ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)((3,3))


<!-- !! processed by numpydoc !! -->

### gaddlemaps.parsers.open_coordinate_file(filename, mode='r')
<!-- !! processed by numpydoc !! -->

* **Return type**

    `CoordinatesParser`
