'''
This module contains the GroFile class and dump_lattice_gro and
extract_lattice_gro functions to avoind cyclic imports.

'''

import warnings
import numpy

from ._itp_parse import (ItpFile, ItpLine, ItpLineAtom, ItpLineBonds,
                         ItpLineMoleculetype, ItpSection)

class GroFile(object):
    """
    Implements a file object created for openning and writting grofiles

    Grofile verifies the gromacs format, and autodetects the correct
    way to read the file. It modifies the methods read and write of a
    common file with methods that take as input lists wth all the info
    of an atom. The readlines and writelines are modified to take
    (return) lists with the info of an atom in lists (as returned by
    the read function).

    This class can also be initiated with an alreadye opened file.

    Note
    ----
    The "+" modifier for the open mode is not allowed right now

    Raises
    ------
    IOError
        If the gro file has not the correct format

    """

    DEFAULT_POSTION_FORMAT = (8, 3)
    DEFAULT_COMMENT = "Gro file genereted with 'Gromacs Tools' python module."
    COORD_START = 20
    NUMBER_FIGURES = 7

    def __init__(self, path, mode="r"):

        super(GroFile, self).__init__()
        
        if hasattr(path, 'read'):
            if path.mode != self._correct_mode(mode):
                raise IOError('Unable to manage an opened file in "{}" mode.'.format(mode))
            mode = path.mode
            self._file = path
        else:
            mode = self._correct_mode(mode)
            self._file = open(path, mode)
        self._comment = None
        self._natoms = None
        self._init_position = None
        self._atomline_bytesize = None
        self._format = {
            "position": None,
            "velocities": None,
        }
        self._box_matrix = numpy.zeros((3, 3))
        self._current_atom = 0

        if "r" in mode:
            self._load_and_verify()

    @property
    def name(self):
        """
        The filename of the open file
        """
        return self._file.name

    @property
    def natoms(self):
        """
        integer: The number of atoms in the system.
        """
        return self._natoms

    @natoms.setter
    def natoms(self, value):
        if "w" in self._file.mode:
            self._natoms = value
        else:
            raise AttributeError("Can't set atributte in read mode")

    @property
    def comment(self):
        """
        string: The comment of the .gro file.
        """
        if self._comment is None:
            return self.DEFAULT_COMMENT
        return self._comment

    @comment.setter
    def comment(self, value):
        if "w" in self._file.mode:
            if value[-1] == '\n':
                value = value[:-1]
            self._comment = value
        else:
            raise AttributeError("Can't set atributte in read mode")

    @property
    def position_format(self):
        """
        tuple : the number of total figures and decimal places for the
        positions
        """
        if self._format["position"] is None:
            return self.DEFAULT_POSTION_FORMAT
        return self._format["position"]

    @position_format.setter
    def position_format(self, value):
        if "w" in self._file.mode:
            self._format["position"] = value
        else:
            raise AttributeError("Can't set atributte in read mode")

    @property
    def box_matrix(self):
        """
        :obj:`numpy.ndarray` (3,3) : the 3 lattice vectors of the box.
        """
        return self._box_matrix

    @box_matrix.setter
    def box_matrix(self, value):
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
            raise IOError("First of a gro line should be not empty")

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

    def seek_atom(self, index):
        """
        Displaces the postion of the 'cursor' to an atom line

        Displaces the postion of the 'cursor' to the begining of the
        line of the 'index' atom, where the first atom index is 0. If
        the index is equal to the number of atoms the begining of the
        box lattice line is found.

        Parameters
        ----------
        index : int
            The index of the atom to found

        """
        self._current_atom = index
        if index > self._natoms:
            error_text = ("Index {} is greater than the number of atoms in the "
                          "grofile").format(index)
            raise IndexError(error_text)

        self._file.seek(self._init_position + index * self._atomline_bytesize)


    def _readline(self):
        return self._file.readline()

    def writelines(self, list_atomlist):
        """
        Writes several lines of atoms

        If there was no content writen in the file it creates the
        header and the number of atoms. If the number of atoms was not
        provided it will kept empty and the number will be writen just
        before closing the file. This posibility is only compatible
        with a number of atoms smaller than 1000000.

        Parameters
        ----------
        atomline : list (list)
            A list of lists with all the info of the atom, just like
            the one returned by realine

        """

        for line in list_atomlist:
            self.writeline(line)

    def writeline(self, atomlist):
        """
        Writes a line of atoms

        If there was no content writen in the file it creates the
        header and the number of atoms. If the number of atoms was not
        provided it will kept empty and the number will be writen just
        before closing the file. This posibility is only compatible
        with a number of atoms smaller than 1000000.

        Parameters
        ----------
        atomlist : list or string
            A list with all the info of the atom, just like the one returned by
            realine. If it is a string, it will be writen directly without
            parsing.

        """
        if self._init_position is None:
            self._setup_write_file(atomlist)
            return

        if isinstance(atomlist, list):
            to_write = self.parse_atomlist(atomlist, format_dict=self._format)
        else:
            to_write = atomlist
        self._file.write(to_write)
        self._file.write("\n")
        self._current_atom += 1

    def _setup_write_file(self, atomlist):
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

    def readline(self, parsed=True):
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
        if self._current_atom >= self._natoms:
            raise StopIteration
        self._current_atom += 1
        return self.parse_atomline(info, self._format)

    def readlines(self):
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

    def next(self):
        """
        Returns next atomline formatted
        """
        return self.readline()

    @classmethod
    def _correct_mode(cls, mode):
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
    def parse_atomlist(cls, atomlist, format_dict=None):
        """
        Parses an atom line and returns its content in a list.

        Parses an atom line and returns its content in a list. In its
        default behaivour the input line is analyzed and the number of
        decimals are guessed. It is also possible to insert the format
        in a format dictionry. This is a dictionary with 2 keys,
        "position" and "velocities". The position key leads to a tuple
        with 2 ints. The first in is the number of characters in the
        float format (C) and the second is the number of decimal
        places in the float (D), i.e. corresponds to the '%C.Df'
        format.

        Parameters
        ----------
        info : list
            A list with the following information:
                + mol_index(integer)
                + resname (string)
                + atomname (string)
                + global_index (integer)
                + x, y, z (floats)
                + vx, vy ,vz (floats) [optional]]

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
        atomline : str
            An string with a line from a grofile
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

        # Wrap the numbers if greater than 99999
        atomlist[0] = atomlist[0] % 99999 + int(atomlist[0] > 99999)
        atomlist[3] = atomlist[3] % 99999 + int(atomlist[3] > 99999)

        # Validate the atomname and resname
        atomlist[1] = cls.validate_string(atomlist[1])
        atomlist[2] = cls.validate_string(atomlist[2])

        format_list = [
            "{:5d}",
            "{:5s}",
            "{:>5s}",
            "{:5d}",
            "{:{figures}.{decimals}f}"*3
        ]
        if velocities:
            format_list.append("{:{figures}.{velocities}f}"*3)
        out = "".join(format_list).format(*atomlist, **float_format_dict)
        return out

    @classmethod
    def parse_atomline(cls, atomline, format_dict=None):
        """
        Parses an atom line and returns its content in a list.

        Parses an atom line and returns its content in a list. In its
        default behaivour the input line is analyzed and the number of
        decimals are guessed. It is also possible to insert the format
        in a format dictionry. This is a dictionary with 2 keys,
        "position" and "velocities". The position key leads to a tuple
        with 2 ints. The first in is the number of characters in the
        float format (C) and the second is the number of decimal
        places in the float (D), i.e. corresponds to the '%C.Df'
        format.

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
                + atomname (string)
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

        info = [
            res_num,
            atomline[5:10].strip(),
            atomline[10:15].strip(),
            atom_num,
            float(atomline[20:20+space]),
            float(atomline[20+space:20+2*space]),
            float(atomline[20+2*space:20+3*space])
        ]

        if format_dict["velocities"]:
            aux_info = [
                float(atomline[20+3*space:20+4*space]),
                float(atomline[20+4*space:20+5*space]),
                float(atomline[20+5*space:20+6*space]),
            ]
            info += aux_info

        return info

    def __next__(self):
        return self.next()

    def __iter__(self):
        while True:
            try:
                yield next(self)
            except StopIteration:
                return

    def close(self):
        """
        Closes the file

        If the filw is in write mode it will write al the remaining information.
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

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @classmethod
    def determine_format(cls, atomline):
        """
        Returns the format of the postion coordinates and wheter or
        not the grofile has velocities

        Parameters
        ----------
        atomline : str
            An atom line from the grofile

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
    def validate_string(string):
        """
        Validates a string to be valid as resname or a as atomname.

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


def _validate_res_atom_numbers(line):
    """
    Validates the residue and number atom in a gro atom_line.
    """
    num_err = "Invalid {} number ({}), check the file format."
    res_num = line[:5]
    atom_num = line[15:20]
    try:
        res_num = int(res_num)
    except ValueError:
        raise IOError(num_err.format('residue', res_num))
    try:
        atom_num = int(atom_num)
    except ValueError:
        raise IOError(num_err.format('atom', atom_num))
    return res_num, atom_num


def extract_lattice_gro(line):
    """
    This function extracts the lattice vectors for pbc from the final
    line of a .gro file.

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


def dump_lattice_gro(vectors):
    """
    This function extracts final line of a .gro file from the lattice vectors
    for pbc.

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

