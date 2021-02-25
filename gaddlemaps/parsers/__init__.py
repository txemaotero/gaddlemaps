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
Simulation Files Parsers
========================

This submodule contains some useful functionalities to read and extract
information from files coming from molecular simulation programs. Up to now,
only files in the format used by GROMACS are contemplated (.gro and .itp
formats).

'''

import warnings
import numpy
import abc
import os.path

from typing import Tuple, List, Union, Dict, Optional, Type

from ._itp_parse import (ItpFile, ItpLine, ItpLineAtom, ItpLineBonds,
                         ItpLineMoleculetype, ItpSection)
from ._top_parsers import read_topology, TopologyParser, TopologyParserManager


GroLine = Union[Tuple[int, str, str, int, float, float, float],
                Tuple[int, str, str, int, float, float, float, float, float, float]]

def open_coordinate_file(filename: str, mode: str="r") -> 'CoordinatesParser':
    name = os.path.basename(filename)
    extension = name.split(".")[-1]

    if extension in ParserManager.parsers:
        return ParserManager.parsers[extension](filename, mode=mode)
    else:
        raise ValueError(f"No parser available for extension {extension}")


class ParserManager:
    parsers: Dict["str", Type["CoordinatesParser"]] = {}

    @classmethod
    def register(cls, parser: Type['CoordinatesParser']):
        if parser.EXTENSIONS:
            for extension in parser.EXTENSIONS:
                cls.parsers[extension] = parser

class ParserRegistered(abc.ABCMeta):
    def __init__(self, name, bases, attrs):
        super().__init__(name, bases, attrs)
        ParserManager.register(self)

    def __new__(metaclass, name, bases, attrs):
        return super().__new__(metaclass, name, bases, attrs)

class CoordinatesParser(metaclass=ParserRegistered):

    EXTENSIONS: Optional[Tuple[str, ...]] = None

    @abc.abstractmethod
    def __init__(self, path: str, mode: str="r"):
        """
        Implements a file object for opening and writing coordinate files

        Grofile verifies the gromacs format, and autodetects the correct
        way to read the file. It modifies the methods read and write of a
        common file with methods that take as input lists with all the info
        of an atom. The readlines and writelines are modified to take
        (return) lists with the info of an atom in lists (as returned by
        the read function).

        This class can also be initiated with an already opened file.

        Note
        ----
        The "+" modifier for the open mode is not allowed right now

        Parameters
        ----------
        path: str
            The path to the file to be opened.
        mode: str
            "r" to open the file in read mode and "w" to open it in write mode.


        Attributes
        ----------
        EXTENSIONS: Optional[Tuple[str, ...]]
            A tuple with all the extensions that the parser is able to process.

        Raises
        ------
        IOError
            If the file has not the correct format
        """

        super().__init__()

    @abc.abstractmethod
    def seek_atom(self, index: int):
        """
        Displaces the position of the 'cursor' to an atom line

        Displaces the position of the 'cursor' to the beginning of the
        line of the 'index' atom, where the first atom index is 0. This function
        warrants that calling the "next" function the information about atom
        "index" will be returned.

        Parameters
        ----------
        index : int
            The index of the atom to found

        """

        return

    @abc.abstractmethod
    def next(self) -> GroLine:
        """
        Returns next atomline formatted
        """
        residue_index = 1
        residue_name = "mol"
        atom_name = "C"
        atom_index = 1
        x_position = 0
        y_position = 0
        z_position = 0

        return (residue_index, residue_name, atom_name, atom_index,
                x_position, y_position, z_position)

    @abc.abstractmethod
    def writeline(self, atomlist: GroLine):
        """
        Writes a line of atom information

        If there was no content written in the file it creates the
        header and the number of atoms. If the number of atoms was not
        provided it will kept empty and the number will be written just
        before closing the file. This possibility is only compatible
        with a number of atoms smaller than 1000000.

        Parameters
        ----------
        atomlist : list or string
            A list with all the info of the atom, just like the one returned by
            readline. If it is a string, it will be written directly without
            parsing.

        """
        return

    @abc.abstractmethod
    def close(self):
        """
        Closes the file
        """
        return

    @property
    @abc.abstractmethod
    def natoms(self) -> int:
        """
        Number of atoms in the file
        """
        return 0

    @property  #type: ignore
    @abc.abstractmethod
    def box_matrix(self) -> numpy.ndarray:
        """
        Return a 3x3 matrix with the 3 lattice vectors
        """
        return numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    @box_matrix.setter  # type: ignore
    @abc.abstractmethod
    def box_matrix(self, new_box: numpy.ndarray):
        """
        This function does no need to work in "r" mode
        """
        return

    @property  #type: ignore
    @abc.abstractmethod
    def comment(self) -> str:
        """
        A comment/name of the system
        """
        return "Generic system"

    @comment.setter  # type: ignore
    @abc.abstractmethod
    def comment(self, new_comment: str):
        """
        This function does not need to work in "r" mode
        """
        return

    def writelines(self, list_atomlist: List[GroLine]):
        """
        Writes several lines of atoms

        Parameters
        ----------
        list_atomlist : list of list of str or int or float
            A list of lists with all the info of the atom, just like
            the one returned by readline

        """

        for line in list_atomlist:
            self.writeline(line)

    def __next__(self):
        return self.next()

    def __iter__(self):
        while True:
            try:
                yield next(self)
            except StopIteration:
                return

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()




class GroFile(CoordinatesParser):
    """
    Implements a file object for opening and writing .gro files

    Grofile verifies the gromacs format, and autodetects the correct
    way to read the file. It modifies the methods read and write of a
    common file with methods that take as input lists with all the info
    of an atom. The readlines and writelines are modified to take
    (return) lists with the info of an atom in lists (as returned by
    the read function).

    This class can also be initiated with an already opened file.

    Note
    ----
    The "+" modifier for the open mode is not allowed right now

    Parameters
    ----------
    path: str
        The path to the file to be opened.
    mode: str
        "r" to open the file in read mode and "w" to open it in write mode.

    Raises
    ------
    IOError
        If the gro file has not the correct format

    """

    EXTENSIONS = ("gro", "GRO")

    DEFAULT_POSTION_FORMAT = (8, 3)
    DEFAULT_COMMENT = "Gro file genereted with 'Gromacs Tools' python module."
    COORD_START = 20

    # Defines the maximum number of atoms in the file 10**NUMBER_FIGURES
    NUMBER_FIGURES = 9

    def __init__(self, path: str, mode: str = "r"):
        mode = self._correct_mode(mode)
        self._file = open(path, mode)
        self._comment: Optional[str] = None
        self._natoms: Optional[int] = None
        self._init_position: Optional[int] = None
        self._atomline_bytesize: Optional[int] = None
        self._format: Dict[str, Optional[Tuple[int, int]]] = {
            "position": None,
            "velocities": None,
        }
        self._box_matrix = numpy.zeros((3, 3))
        self._current_atom = 0

        if "r" in mode:
            self._load_and_verify()

    @property
    def name(self) -> str:
        """
        str: The filename of the open file
        """
        return self._file.name

    @property
    def natoms(self) -> int:
        """
        integer: The number of atoms in the system.
        """
        if self._natoms is None:
            raise ValueError("The system was not fully initialized")
        return self._natoms

    @natoms.setter
    def natoms(self, value: int):
        if "w" in self._file.mode:
            self._natoms = value
        else:
            raise AttributeError("Can't set atributte in read mode")

    @property
    def comment(self) -> str:
        """
        string: The comment of the .gro file.
        """
        if self._comment is None:
            return self.DEFAULT_COMMENT
        return self._comment

    @comment.setter
    def comment(self, value: str):
        if "w" in self._file.mode:
            if value[-1] == '\n':
                value = value[:-1]
            self._comment = value
        else:
            raise AttributeError("Can't set atributte in read mode")

    @property
    def position_format(self) -> Tuple[int, int]:
        """
        tuple of int : the number of total figures and decimal places for the
            positions.
        """
        if self._format["position"] is None:
            return self.DEFAULT_POSTION_FORMAT
        return self._format["position"]

    @position_format.setter
    def position_format(self, value: Tuple[int, int]):
        if "w" in self._file.mode:
            self._format["position"] = value
        else:
            raise AttributeError("Can't set atributte in read mode")

    @property
    def box_matrix(self) -> numpy.ndarray:
        """
        numpy.ndarray(3,3): the 3 lattice vectors of the box.
        """
        return self._box_matrix

    @box_matrix.setter
    def box_matrix(self, value: numpy.ndarray):
        if "r" in self._file.mode:
            raise AttributeError("Can't set atributte in read mode")
        value = numpy.array(value)

        shape = value.shape
        if shape == (3,):
            value = numpy.diag(value)
        elif shape == (3, 3):
            pass
        else:
            text = ("The input array have {} dimension. Input arrays must be a "
                    "3x3 array or a vector with 3 elements.").format(shape)
            raise ValueError(text)

        self._box_matrix = value

    def _load_and_verify(self):
        """
        Verify if the file can be a grofile and determine the info in it
        """

        self._comment = self._readline()

        if not self._comment:
            raise IOError("First line of a gro line must not be empty")

        try:
            self._natoms = int(self._readline())
        except ValueError:
            raise IOError(("Second line of a gro file should contain only "
                           "the number of atoms and blank spaces"))

        # Save this position to ve able to return later
        self._init_position = self._file.tell()

        # read the first line of atoms
        first_line = self._readline()
        self._format = self.determine_format(first_line)
        self._atomline_bytesize = self._file.tell() - self._init_position
        self._load_box_matrix()
        # Return to the beggining of the atoms
        self.seek_atom(0)

    def _load_box_matrix(self):
        """Seek after the last atom and loads the box info"""
        self.seek_atom(self._natoms)
        line = self._readline()
        if not line:
            error_text = ("The number of atoms at the top of the file does not "
                          "correspond to the actual number of atoms or the "
                          "simulation box line is missing from the file")
            raise IOError(error_text)
        # Try to load the simulation box
        try:
            self._box_matrix = extract_lattice_gro(line)
        except ValueError:
            error_text = ("The number of atoms at the top of the file does not "
                          "correspond to the actual number of atoms or the "
                          "simulation box line has the wrong fomat")
            raise IOError(error_text)

    def seek_atom(self, index: int):
        """
        Displaces the position of the 'cursor' to an atom line

        Displaces the position of the 'cursor' to the beginning of the
        line of the 'index' atom, where the first atom index is 0. If
        the index is equal to the number of atoms the beginning of the
        box lattice line is found.

        Parameters
        ----------
        index : int
            The index of the atom to found

        """
        self._current_atom = index
        if index > self.natoms:
            error_text = ("Index {} is greater than the number of atoms in the "
                          "grofile").format(index)
            raise IndexError(error_text)

        if self._init_position is None or self._atomline_bytesize is None:
            raise ValueError("Error in initalizaiton")

        self._file.seek(self._init_position + index * self._atomline_bytesize)

    def _readline(self) -> str:
        return self._file.readline()

    def writelines(self, list_atomlist: List[GroLine]):
        """
        Writes several lines of atoms

        If there was no content written in the file it creates the
        header and the number of atoms. If the number of atoms was not
        provided it will kept empty and the number will be written just
        before closing the file. This possibility is only compatible
        with a number of atoms smaller than 1000000.

        Parameters
        ----------
        list_atomlist : list of list of str or int or float
            A list of lists with all the info of the atom, just like
            the one returned by readline

        """

        for line in list_atomlist:
            self.writeline(line)

    def writeline(self, atomlist: GroLine):
        """
        Writes a line of atom information

        If there was no content written in the file it creates the
        header and the number of atoms. If the number of atoms was not
        provided it will kept empty and the number will be written just
        before closing the file. This possibility is only compatible
        with a number of atoms smaller than 1000000.

        Parameters
        ----------
        atomlist : list or string
            A list with all the info of the atom, just like the one returned by
            readline. If it is a string, it will be written directly without
            parsing.

        """
        if self._init_position is None:
            self._setup_write_file(atomlist)
            return

        if isinstance(atomlist, (list, tuple)):
            to_write = self.parse_atomlist(atomlist, format_dict=self._format)
        else:
            to_write = atomlist
        self._file.write(to_write)
        self._file.write("\n")
        self._current_atom += 1

    def _setup_write_file(self, atomlist: Union[str, GroLine]):
        if isinstance(atomlist, str):
            atomlist = self.parse_atomline(atomlist)
        # setupt the format_dict
        if self._format["position"] is None:
            self._format["position"] = self.DEFAULT_POSTION_FORMAT
            self._format["velocities"] = len(atomlist) == 10

        # There is no content writen to the file
        comment = self.comment
        self._file.write(comment)
        if comment[-1] != "\n":
            self._file.write("\n")

        # Write the number of atoms if possible
        if self._natoms is None:
            self._file.write(" "*self.NUMBER_FIGURES+"\n")
        else:
            self._file.write("{:d}\n".format(self._natoms))
        self._init_position = self._file.tell()
        self.writeline(atomlist)
        self._atomline_bytesize = self._file.tell() - self._init_position

    def readline(self, parsed: bool=True) -> Union[GroLine, str]:
        """
        Returns the next line of the gro file.

        If parsed=True (default) the data on the atom line is parsed
        and returned in a list. If false the line is returned as an
        string (just like if it were opened as a common file.)

        Parameters
        ----------
        parsed : bool (optional)
            Wether the data of the line should be parsed or not.
            Default: True

        Returns
        -------
        info : str or list
            Returns the info of the line, i alist if it was parsed or
            in a string if not

        """
        info = self._readline()
        if not parsed:
            return info
        if self._current_atom >= self.natoms:
            raise StopIteration
        self._current_atom += 1
        return self.parse_atomline(info, self._format)

    def readlines(self) -> List[GroLine]:
        """
        Returns all the lines of a grofile parsed in lists.

        Returns
        -------
        info : list(list)
            Returns a list with the info of all the atom lines in lists

        """
        # self.seek_atom(0)
        info = [i for i in self]
        return info

    def next(self) -> GroLine:
        """
        Returns next atomline formatted
        """
        return self.readline()  #type: ignore

    @classmethod
    def _correct_mode(cls, mode: str) -> str:
        """
        Checks if the '+' modifier is set and removes it because it is
        not compatible yet
        """
        if "+" in mode:
            mode = mode.replace("+", "")
            mode = mode.replace("w", "r")
            warn_text = ("'+' modifier for open mode is not implemented, the"
                         " file will be opened in {} mode").format(mode)
            warnings.warn(warn_text, RuntimeWarning)
        return mode

    @classmethod
    def parse_atomlist(cls, atomlist: GroLine,
                       format_dict: Dict = None) -> str:
        """
        Convert a list of atom info to string with the appropriate format

        Parses a list with the atom information and returns its content in a
        string with the correct format to write in a .gro file. In its default
        behaviour the input line is analyzed and the number of decimals are
        guessed. It is also possible to insert the format in a format
        dictionary. This is a dictionary with 2 keys, "position" and
        "velocities". The position key leads to a tuple with 2 ints. The
        first is the number of characters in the float format (C) and the
        second is the number of decimal places in the float (D), i.e.
        corresponds to the '%C.Df' format.

        Parameters
        ----------
        atomlist: list
            A list with the following information:
                + mol_index(integer)
                + resname (string)
                + name (string)
                + global_index (integer)
                + x, y, z (floats)
                + vx, vy ,vz (floats) [optional]]

        format_dict: dict (optional)
            None or dictionary with 2 keys, "position" and "velocities". The
            position key leads to a tuple with 2 ints. THe first in is
            the number of characters in the float format (C) and the
            second is the number of decimal places in the float (D),
            i.e. corresponds to the '%C.Df' format. If None the format
            will be automatically guessed.
            Default None.

        Returns
        -------
        atomline : str
            An string with a line from a .gro file
        """

        if format_dict is None:
            float_format = cls.DEFAULT_POSTION_FORMAT
        else:
            float_format = format_dict["position"]

        float_format_dict = {
            "figures": float_format[0],
            "decimals": float_format[1]
        }

        if len(atomlist) == 10:
            velocities = True
            float_format_dict["velocities"] = float_format_dict["decimals"]+1
        elif len(atomlist) == 7:
            velocities = False
        else:
            length = len(atomlist)
            text = "Atomlist must have len 7 or 10, found {}".format(length)
            raise ValueError(text)

        # If an external format was set check if the velocities
        # requirement is met

        if format_dict is not None:
            if velocities == format_dict["velocities"]:
                pass
            elif velocities:
                text = ("Error while parsing {}. The atom information has "
                        "velocities, but the first atom of the GroFile has not."
                        " All the atoms must have the number of fields"
                        " to write.").format(atomlist)
                raise IOError(text)
            else:
                text = ("Error while parsing {}. The atom information has not"
                        "velocities, but the first atom of the GroFile has."
                        " All the atoms must have the number of fields"
                        " to write.").format(atomlist)
                raise IOError(text)

        atominfo = list(atomlist)
        # Wrap the numbers if greater than 99999
        atominfo[0] = atomlist[0] % 99999 + int(atomlist[0] > 99999)
        atominfo[3] = atomlist[3] % 99999 + int(atomlist[3] > 99999)

        # Validate the name and resname
        atominfo[1] = cls.validate_string(atomlist[1])
        atominfo[2] = cls.validate_string(atomlist[2])

        format_list = [
            "{:5d}",
            "{:5s}",
            "{:>5s}",
            "{:5d}",
            "{:{figures}.{decimals}f}"*3
        ]
        if velocities:
            format_list.append("{:{figures}.{velocities}f}"*3)
        out = "".join(format_list).format(*atominfo, **float_format_dict)
        return out

    @classmethod
    def parse_atomline(cls, atomline: str,
                       format_dict: Optional[Dict] = None) -> GroLine:
        """
        Parses an atom line and returns its content in a list.

        Parses an atom line and returns its content in a list. In its default
        behaviour the input line is analyzed and the number of decimals are
        guessed. It is also possible to insert the format in a format
        dictionary. This is a dictionary with 2 keys, "position" and
        "velocities". The position key leads to a tuple with 2 ints. The
        first in is the number of characters in the float format (C) and the
        second is the number of decimal places in the float (D), i.e.
        corresponds to the '%C.Df' format.

        Parameters
        ----------
        atomline : str
            An string with a line from a grofile
        format_dict : dict_format (optional)
            None or dictionary with 2 keys, "position" and "velocities". The
            position key leads to a tuple with 2 ints. THe first in is
            the number of characters in the float format (C) and the
            second is the number of decimal places in the float (D),
            i.e. corresponds to the '%C.Df' format. If None the format
            will be automatically guessed.
            Default None.

        Returns
        -------
        info : list
            A list with the following information:
                + mol_index(integer)
                + resname (string)
                + name (string)
                + global_index (integer)
                + x, y, z (floats)
                + vx, vy ,vz (floats) [optional]

        """

        if format_dict is None:
            format_dict = cls.determine_format(atomline)

        if atomline[-1] == "\n":
            atomline = atomline[:-1]

        length = len(atomline)
        space = format_dict["position"][0]
        expected_length = 20+space*3*(1+format_dict["velocities"])

        if length != expected_length:
            texterr = ("Found an atom line with {} characters, when {} were "
                       "expected. Check your file.").format(length,
                                                            expected_length)
            raise IOError(texterr)

        res_num, atom_num = _validate_res_atom_numbers(atomline)

        info: GroLine = (
            res_num,
            atomline[5:10].strip(),
            atomline[10:15].strip(),
            atom_num,
            float(atomline[20:20+space]),
            float(atomline[20+space:20+2*space]),
            float(atomline[20+2*space:20+3*space])
        )

        if format_dict["velocities"]:
            aux_info: Tuple[float, float, float] = (
                float(atomline[20+3*space:20+4*space]),
                float(atomline[20+4*space:20+5*space]),
                float(atomline[20+5*space:20+6*space]),
            )
            info += aux_info  #type: ignore

        return info

    def close(self):
        """
        Closes the file

        If the file is in write mode it will write al the remaining
        information.
        """
        if "w" in self._file.mode:
            self._write_closing_info()

        self._file.close()

    def _write_closing_info(self):
        """
        Writes alll the necesary info when the file is in write mode
        """
        if self._natoms is None and self._current_atom == 0:
            warnings.warn("Closing an empty file")
            return
        # If no number of atoms was given at the begining fill the value
        if self._natoms is None:
            self._file.seek(self._init_position-1-self.NUMBER_FIGURES)
            self._natoms = self._current_atom
            self._file.write("{:{figures}d}\n".format(self._current_atom,
                                                      figures=self.NUMBER_FIGURES))

        # If there was a number, check that it is the correct one
        else:
            if self._natoms != self._current_atom:
                text = ("{} Atoms were writen but {} were expected.")
                text = text.format(self._current_atom, self._natoms)
                raise IOError(text)


        self.seek_atom(self._natoms)
        self._file.write(dump_lattice_gro(self._box_matrix))
        self._file.write("\n")



    @classmethod
    def determine_format(cls, atomline: str) -> Dict:
        """
        Returns the format of the postion coordinates and whether or
        not the .gro file has velocities

        Parameters
        ----------
        atomline : str
            An atom line from the .gro file

        Returns
        -------
        form : dict
            Dictionary with 2 keys, "position" and "velocities". The
            position key leads to a tuple with 2 ints. THe first in is
            the number of characters in the float format (C) and the
            second is the number of decimal places in the float (D),
            i.e. corresponds to the '%C.Df' format.
        """

        # remove \n at the end if present
        if atomline[-1] == "\n":
            atomline = atomline[:-1]
        size = len(atomline)

        # Check that thereb are not more "\n" if the atomline
        if atomline.count("\n"):
            raise ValueError("The input must correspond to only 1 line")

        # count the number of "." fater the coordinate start
        ndots = atomline[cls.COORD_START:].count(".")

        if ndots == 3:
            velocities = False
        elif ndots == 6:
            velocities = True
        else:
            errortext = ("Found {} '.' after the 20 first characters. This is "
                         "not compatible with the gro format. "
                         "Check the input file").format(ndots)
            raise IOError(errortext)

        nfigures = (size - cls.COORD_START) // ndots

        if size != (cls.COORD_START + ndots * nfigures):
            raise IOError(("Some of the coordinates/velocities values has "
                           "different number of characters. "
                           "This is not allowed"))
        ndecimals = nfigures - 5

        form = {
            "position": (nfigures, ndecimals),
            "velocities" : velocities,
        }

        return form

    @staticmethod
    def validate_string(string: str) -> str:
        """
        Validates a string to be valid as resname or as name.

        If the input name has more than 5 characters, new name is returned
        cutting the input one.

        Parameters
        ----------
        string : str
            The string to validate.

        Returns
        -------
        new_string : str
            The properly formatted string.

        """
        if len(string) > 5:
            old_string = string
            string = string[:5]
            warn_text = ('The input resname has more than 5 character ({})'
                         ', this is modified to be {}.').format(old_string,
                                                                string)
            warnings.warn(warn_text)
        return string


def _validate_res_atom_numbers(line: str) -> Tuple[int, int]:
    """
    Validates the residue and number atom in a gro atom_line.
    """
    num_err = "Invalid {} number ({}), check the file format."
    res_num_str = line[:5]
    atom_num_str = line[15:20]
    try:
        res_num = int(res_num_str)
    except ValueError:
        raise IOError(num_err.format('residue', res_num_str))
    try:
        atom_num = int(atom_num_str)
    except ValueError:
        raise IOError(num_err.format('atom', atom_num_str))
    return res_num, atom_num


def extract_lattice_gro(line: str) -> numpy.ndarray:
    """
    Extracts the lattice vectors for pbc from the final line of a .gro file.

    Parameters
    ----------
    line : string
        The final line of a .gro file.

    Returns
    -------
    vectors : numpy.ndarray((3,3))
        An array with a lattice vector in each row.

    """

    vectors = numpy.zeros(9)
    index = (0, 4, 8, 1, 2, 3, 5, 6, 7)
    nums = (float(i) for i in line.split())
    for i, num in zip(index, nums):
        vectors[i] = num
    return vectors.reshape((3, 3))


def dump_lattice_gro(vectors: numpy.ndarray) -> str:
    """
    Extracts final line of a .gro file from the lattice vectors for pbc.

    Parameters
    ----------
    vectors : numpy.ndarray((3,3))
        An array with a lattice vector in each row.

    Returns
    -------
    line : string
        The final line of a .gro file.

    """

    # Reorder the vectors comps
    vectors = vectors.ravel()
    index = (0, 4, 8, 1, 2, 3, 5, 6, 7)
    new_vector = numpy.zeros_like(vectors)
    for i, val in enumerate(index):
        new_vector[i] = vectors[val]
    if numpy.any(new_vector[3:]):
        lim = 9
    else:
        lim = 3
    return ' '.join(['{:9.5f}'.format(i) for i in new_vector[:lim]])


