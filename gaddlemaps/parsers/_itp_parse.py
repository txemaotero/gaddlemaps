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
This module contains useful features to parse itp files.
'''


import re
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Tuple
from warnings import warn


class ItpFile(OrderedDict):
    """
    A class to handle itp files based on dictionary.

    The keys are the titles of the found sections and the values are ItpSection
    objects. There is the 'header' key that contains a list of the lines
    before any section.

    Parameters
    ----------
    fitp : string
        The itp file to handle.

    """
    def __init__(self, fitp: str):
        super(ItpFile, self).__init__()
        self.fitp = fitp

        self._parse_file()

    def _parse_file(self):
        sec = None
        self['header'] = []
        _file = open(self.fitp, encoding='utf-8')

        for line in _file:
            # If not section line
            if not re.match(r'\[.*\]', line.strip()):
                if sec is None:
                    self['header'].append(line)
                else:
                    self[sec].append(line)
            else:
                sec = re.findall(r'\[(.*)\]', line)[0].strip()
                self[sec] = ItpSection(sec, [])
        _file.close()

    def __str__(self):
        msg = 'Itp file with the following sections:\n'
        msg += '\n'.join(['{}'.format(sec) for sec in sorted(self.keys())])
        return msg

    def __repr__(self):
        return '"{}" itp file.\n'.format(self.fitp)

    def copy(self) -> 'ItpFile':
        """
        Returns a copy of the object that is not self.

        Returns
        -------
        itp_file : ItpFile
            The self copy.
        """
        return ItpFile(self.fitp)

    def write(self, fout: str):
        """
        Writes the itp content in fout.

        Parameters
        ----------
        fout : str
            The path with the output itp.

        """
        fopen = open(fout, 'w')
        for header, section in self.items():
            if header == 'header':
                for line in section:
                    fopen.write(line)
            else:
                fopen.write('{}\n'.format(section))


class ItpSection(list):
    """
    A class to wrap .itp sections based on list.

    The object only contains lines with content but the whole section is
    accessible through self.lines.

    Parameters
    ----------
    section_name : string
        The header of the section.
    lines : list of string
        A list with the lines of the section.

    Raises
    ------
    IOError
        If one of the input lines corresponds to a section header.

    """
    def __init__(self, section_name: str, lines: List[str]):
        super(ItpSection, self).__init__()
        self._section_name = section_name
        self._lines: List['ItpLine'] = []
        for line in lines:
            self.append(line)

    @staticmethod
    def parse_line(line: str, sec_name: str) -> 'ItpLine':
        """
        Returns the corresponding ItpLine object depending on the section name.

        Parameters
        ----------
        line : string
            The line from itp file to parse.
        sec_name : string
            The itp section name.

        Returns
        -------
        itp_line : ItpLine
            An instance of ItpLine, different for each section: ItpLineAtom,
            ItpLineBonds, ItpLineMoleculetype, ItpLine.

        """
        if sec_name == 'atoms':
            return ItpLineAtom(line)
        elif sec_name in ['bonds', 'constraints', 'pairs']:
            return ItpLineBonds(line)
        elif sec_name in 'moleculetype':
            return ItpLineMoleculetype(line)
        return ItpLine(line)

    def append(self, new_line: str):
        parse = self.parse_line(new_line, self._section_name)
        if parse.content:
            super(ItpSection, self).append(parse)
        self._lines.append(parse)

    def __str__(self):
        msg = '[ {} ]\n'.format(self.section_name)
        msg += ''.join((str(line) for line in self._lines))
        return msg

    def __repr__(self):
        return 'Itp "{}" section.'.format(self.section_name)

    @property
    def section_name(self) -> str:
        """
        string : The header of the section
        """
        return self._section_name

    @property
    def lines(self) -> List['ItpLine']:
        """
        list of ItpLine : A list with every line of the section.
        """
        return self._lines[:]


class ItpLine:
    """
    A class to parse content of lines from .itp file.

    This class is used as parent for more specific types of lines.

    Parameters
    ----------
    line : string
        The itp line to parse.

    Attributes
    ----------
    content : string
        The important information of the line, i.e. atom info.
    comment : string
        The comment in the line. If there is no comment it is set to empty
        str.

    Raises
    ------
    IOError
        If the input line corresponds to a section header.

    """
    def __init__(self, line: str):
        super(ItpLine, self).__init__()

        self._content, self._comment = self.parse_itp_line(line)

    def __str__(self):
        return self.line

    __repr__ = __str__

    @property
    def line(self) -> str:
        """
        string : the complete content of the itp line.
        """
        if self.comment:
            if self.comment.startswith('#'):
                return self._comment
            return '; '.join((self._content, self._comment))
        return self._content

    @property
    def content(self) -> str:
        """
        string : The relevant information of the line (no comments).
        """
        return self._content.strip()

    @content.setter
    def content(self, new_content: str):
        warn(('Changing the content of an itp line can lead to compilation '
              'errors'))
        self._content = new_content

    @property
    def comment(self) -> str:
        """
        string : The comment of the line.
        """
        return self._comment.strip()

    @comment.setter
    def comment(self, new_comment: str):
        self._comment = new_comment

    @classmethod
    def parse_itp_line(cls, line: str) -> Tuple[str, str]:
        """
        Extract the content and the comment from an .itp line.

        If line starts with "#" it will be treated as a comment line.

        Parameters
        ----------
        line : string
            The itp line to parse.

        Returns
        -------
        content : string
            The important information of the line, i.e. atom info.
        comment : string
            The comment in the line. If there is no comment it is set to empty
            str.

        Raises
        ------
        IOError
            If the input line corresponds to a section header.

        """
        if not line.strip():
            return '', ''

        if re.match(r'\[.*\]', line):
            raise IOError(('The input line is a section header, can not be '
                           'parsed: {}').format(line))

        if line.startswith('#'):
            return '', line[:]

        if line.startswith(';'):
            return '', line[1:]

        if ';' in line[:-1]:
            spl = line.split(';')
            return spl[0], ';'.join(spl[1:])

        if line[-1] == ';':
            return line[:-1], ''

        # The case with no comments
        return line, ''


class ItpLineAtom(ItpLine):
    """
    A class to parse atom content line from .itp file, based on ItpLine.

    Parameters
    ----------
    line : string
        The itp line to parse.

    """
    def __init__(self, line: str):
        super(ItpLineAtom, self).__init__(line)
        self._fields: Dict[str, Any] = {}
        if self.content:
            self._init_fields()

    def _init_fields(self):
        parsed_content = self.content.split()
        self._fields['nr'] = int(parsed_content[0])
        self._fields['type'] = parsed_content[1]
        self._fields['resid'] = int(parsed_content[2])
        self._fields['resname'] = parsed_content[3]
        self._fields['name'] = parsed_content[4]
        self._fields['cgnr'] = int(parsed_content[5])
        try:
            self._fields['charge'] = float(parsed_content[6])
        except IndexError:
            self._fields['charge'] = None
        try:
            self._fields['mass'] = float(parsed_content[7])
        except IndexError:
            self._fields['mass'] = None
        if len(parsed_content) > 8:
            for i, field in enumerate(parsed_content[8:]):
                self._fields['extra{}'.format(i)] = field

    @property
    def number(self) -> int:
        """
        int : The atom number in the itp.
        """
        return self._fields['nr']

    @property
    def type(self) -> str:
        """
        string : The atom type in the itp.
        """
        return self._fields['type']

    @property
    def resid(self) -> int:
        """
        int : The atom residue number in the itp.
        """
        return self._fields['resid']

    @property
    def resname(self) -> str:
        """
        string : The atom residue name in the itp.
        """
        return self._fields['resname']

    @property
    def name(self) -> str:
        """
        string : The atom name in the itp.
        """
        return self._fields['name']

    @property
    def cgnr(self) -> int:
        """
        int : The atom charge group number in the itp.
        """
        return self._fields['cgnr']

    @property
    def charge(self) -> float:
        """
        float : The atom charge in the itp.
        """
        return self._fields['charge']

    @charge.setter
    def charge(self, val: float):
        self._fields["charge"] = val

    @property
    def mass(self) -> float:
        """
        float : The atom mass in the itp.
        """
        mass = self._fields['mass']
        if mass is None:
            raise AttributeError('The mass for this atom is not defined')
        return self._fields['mass']

    @mass.setter
    def mass(self, val: float):
        self._fields["mass"] = val

    @property
    def pair_interaction(self):
        """
        The object storing the pair interaction
        """
        if "pair_interaction" not in self._fields:
            text = ("The pair interaction is not defined. Use the "
                    "ForceFieldData object to generate it.")
            raise AttributeError(text)
        return self._fields["pair_interaction"]

    @pair_interaction.setter
    def pair_interaction(self, val):
        self._fields["pair_interaction"] = val

    @property
    def aux_name(self) -> str:
        """
        An auxiliar name derived from the forcefield
        """
        if "aux_name" not in self._fields:
            text = ("The auxiliar name is not defined. Use the "
                    "ForceFieldData object to generate it.")
            raise AttributeError(text)
        return self._fields["aux_name"]

    @aux_name.setter
    def aux_name(self, val: str):
        self._fields["aux_name"] = val

    @property
    def parsed_line(self) -> List:
        """
        List of int or str: Return a list with the found fields in the atom
            line.
        """
        lis = [
            self._fields['nr'],
            self._fields['type'],
            self._fields['resid'],
            self._fields['resname'],
            self._fields['name'],
            self._fields['cgnr'],
            self._fields['charge'],
        ]
        if self._fields['mass'] is not None:
            lis.append(self._fields['mass'])
        return lis

    @classmethod
    def read_itp_atom(cls, line: str) -> List[Any]:
        return cls(line).parsed_line


class ItpLineBonds(ItpLine):
    """
    A class to parse atom bonds or pairs line from .itp file, based on ItpLine.

    Parameters
    ----------
    line : string
        The itp line to parse.

    """
    def __init__(self, line: str):
        super(ItpLineBonds, self).__init__(line)
        self._fields: Dict[str, Any] = {}
        if self.content:
            self._init_fields()

    def _init_fields(self):
        fields = self.content.split()
        n_fieds = len(fields)
        if n_fieds < 2:
            raise IOError('Wrong format for a bonds line.')
        self._fields['ai'] = int(fields[0])
        self._fields['aj'] = int(fields[1])
        try:
            self._fields['funct'] = int(fields[2])
        except IndexError:
            self._fields['funct'] = 1
        self._fields['const'] = aux = []
        for i_field in range(3, n_fieds):
            try:
                cons0 = float(fields[i_field])
            except ValueError:
                cons0 = fields[i_field]
            aux.append(cons0)

    @property
    def atom_from(self) -> int:
        """
        int : The index of the first bonded atom (ai).
        """
        return self._fields['ai']

    @property
    def atom_to(self) -> int:
        """
        int : The index of the second bonded atom (aj).
        """
        return self._fields['aj']

    @property
    def funct(self) -> int:
        """
        int : The function of the bond.
        """
        return self._fields['funct']

    def __getattr__(self, const):
        if const.startswith('const'):
            try:
                num = int(const[5:])
                val = self._fields['const'][num]
            except (ValueError, IndexError):
                pass
            else:
                return val
        raise AttributeError(('{} has not attribute {}'
                              '').format(self.__class__.__name__, const))


class ItpLineMoleculetype(ItpLine):
    """
    A class to parse the atom moleculetype from .itp file, based on ItpLine.

    Parameters
    ----------
    line : string
        The itp line to parse.

    """
    def __init__(self, line: str):
        super(ItpLineMoleculetype, self).__init__(line)
        if self.content:
            parse = line.split()
            self.name = parse[0]
            self.nrexcl = int(parse[1])
        else:
            self._name: Optional[str] = None
            self._nrexcl: Optional[int] = None

    @property
    def name(self) -> Optional[str]:
        """
        string : The molecule name.
        """
        return self._name

    @name.setter
    def name(self, new_name: str):
        if not isinstance(new_name, str):
            raise TypeError('The molecule name has to be a string.')
        if ' ' in new_name:
            raise ValueError('The molecule name must not contain spaces.')
        self._name = new_name

    @property
    def nrexcl(self) -> Optional[int]:
        """
        int : The number of neighbour atoms to account for the interactions.
        """
        return self._nrexcl

    @nrexcl.setter
    def nrexcl(self, new_val: int):
        if not isinstance(new_val, int):
            raise TypeError('nrexcl has to be an integer.')
        if new_val < 1:
            raise ValueError('nrexcl hast to be equal to or greater than one.')
        self._nrexcl = new_val
