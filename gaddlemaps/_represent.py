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
from . import Alignment, Manager, Restriction
from .components import Residue
from .parsers import GroFile
from ._manager import Restrictions
from typing import List, Tuple, Dict, Optional
from ipywidgets import Button, HBox, VBox, Layout, Label, Box, Tab, Accordion, Widget
from os.path import join, isdir
from os import mkdir




class Cycle:
    """
    A simply generator similar to itertools.cycle, but with the ability to reset
    it to the starting configuration
    """
    def __init__(self, values):
        self._values = values
        self._index = 0

    def __iter__(self):
        while True:
            yield self._values[self._index]
            self._index += 1
            self._index %= len(self._values)

    def __next__(self):
        val = self._values[self._index]
        self._index += 1
        self._index %= len(self._values)
        return val

    def reset(self):
        self._index = 0

def nglview_struct(molecule: Residue):
    """
    Creates an object that acts as a interface betweeen gaddle maps residues and
    nglview. The class is stored inside a funciton to isolate the dependecy.

    Parameters
    ----------
    molecule: Residue
        The molecule that will be represented in nglview

    Returns
    -------
    structure: nglview.Structure
        A structure that can be represented in nglview

    """
    import nglview
    class NglResidue(nglview.Structure):
        def __init__(self, molecule: Residue, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self._molecule = molecule
                self.ext = "gro"

        def get_structure_string(self):
            header = f"Generated for using with GADDLE Maps\n {len(self._molecule)}\n"
            footer = f"0.0 0.0 0.0\n"

            data = "\n".join([GroFile.parse_atomlist(atom.gro_line()) for atom in self._molecule])

            return header + data + footer

    return NglResidue(molecule)

def create_widget_restrictions(mol_low_res: Residue, mol_high_res: Residue,
                               restrict: Restriction = None)-> Tuple[Box, Restriction]:
    """
    Creates a jupyter widget that eases the creation of restraints between 2
    resolution of molecules.

    Parameters
    ----------
    mol_low_res: Residue
        The molecule in the starting reolution
    mol_high_res: Residue
        The moleucle in the end resolution
    restrric: Optional[List[Tuple[int, int]]]:
        If not None is a list where the new restrictions will be appended.

    Returns
    -------
    box: ipywidgets.Box
        The box taht contains all the necessary widgets ready to be displayed

    restrictions: List[Tuple[int, int]]
        The list where the generated restrictions will be stored. If restric was
        not None, this is the same object than the input
    """
    import nglview

    if restrict is not None:
        restrictions = restrict
    else:
        restrictions = []



    text = Label(value="Low Resolution")
    view = nglview.NGLWidget(nglview_struct(mol_low_res))

    text2 = Label(value="High Resolution")
    view2 = nglview.NGLWidget(nglview_struct(mol_high_res))

    button = Button(description="Add restriction")
    button2 = Button(description="Delete restriction")

    restrictions_text = Label(value="Restrictions =")
    restrictions_val = Label(value="")

    colors = Cycle(["red", "blue", "green", "white", "grey", "#ff00ff"])

    def update_restrictions(reset: bool=False):
        if reset:
            view.clear_representations()
            view2.clear_representations()
            colors.reset()

            for restriction in restrictions:
                color = next(colors)
                view.add_representation("licorice", radius=1, color=color,
                                        selection=[restriction[0]], opacity=0.75)
                view2.add_representation("licorice", radius=1, color=color,
                                         selection=[restriction[1]], opacity=0.75)


        view.add_representation('licorice', selection="all")
        view2.add_representation('licorice', selection="all")
        restrictions_val.value = ", ".join([str(i) for i in restrictions])

    def add_restriction(*args):
        color = next(colors)

        restrictions.append((view.picked["atom1"]["index"], view2.picked["atom1"]["index"]))

        view.add_representation("licorice", radius=1, color=color, selection=[restrictions[-1][0]], opacity=0.75)
        view2.add_representation("licorice", radius=1, color=color, selection=[restrictions[-1][1]], opacity=0.75)

        update_restrictions()

    def delete_restriction(*args):
        if restrictions:
            restrictions.pop()
        update_restrictions(reset=True)

    button.on_click(add_restriction)
    button2.on_click(delete_restriction)
    update_restrictions(reset=True)
    box = VBox([text, view, text2, view2, HBox([button, button2]), HBox([restrictions_text, restrictions_val])])
    return box, restrictions


def create_interactive_restriction(correspondence: Dict[str, Alignment]) -> Tuple[Dict[str, Box],
                                                                                  Restrictions]:
    """
    Initialices all the Boxes with the widgets and the dictionary with the
    restrictions for a given manager.

    Parameters
    ----------
    correspondence: Dict[str, Alignment]
        A dictionary that takes as key a molecule name and its aligmnet as value.

    Returns:
    --------
    boxes: Dict[str, ipywidgets.Box]
        A dictionary that for each species returns the box with the widget
        initialized for creating the restrictions.
    restrictions: Dict[str, List[Tuple[int, int]]]
        The dictionary with the restrictions that will be generated by the
        widget for each specie. This will be initially empty and it will be
        filled as the widget is used.
    """

    restrictions = {}
    boxes = {}

    for specie in correspondence:
        try:
            box, restrict = correspondence[specie].interactive_restrictions()
        except OSError:
            continue

        boxes[specie] = box
        restrictions[specie] = restrict

    return boxes, restrictions

def interactive_restrictions(correspondence: Dict[str, Alignment],
                             style:int=None) -> Tuple[Widget, Restrictions]:
    """
    Creates the widget to generate the restrictions of all the species in the
    alignment. It generates the final representation fo the widget.

    Parameters
    ----------

    manager: Manager
        The object that manages all the alignment process
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

    if style is None:
        style = 2

    if style == 0:
        representation = Tab
    elif style == 1:
        representation = Accordion
    elif style == 2:
        representation = VBox
    else:
        representation = VBox

    boxes, restrictions = create_interactive_restriction(correspondence)

    restriction_widget = representation()

    for index, specie in enumerate(boxes):
        child = boxes[specie]
        if style in [2, ]:
            child.children = (Label(value=r"$\textbf{"+ specie + r"}$"), ) + child.children
        restriction_widget.children += (boxes[specie],)
        if style not in [2, ]:
            restriction_widget.set_title(index, specie)

    return restriction_widget, restrictions


def compare_molecules(molecule_low_res: Residue, molecule_high_res: Residue,
                        radius: float=None):
    """
    Creates a visualtization of two molecules to compare their similarities. The
    first molecule is represented with low opacity and with big chunky atoms,
    while the second one is represented in the standard fashion.

    Parameters
    ----------
    molecule_low_res: Residue
        The molecule in the low resolution
    molecule_high_res: Residue
        The molecule in the high resolution
    radius: Optional[float]
        The radius of the atoms in the low resolution representation. Default
        2.5

    Returns
    -------
    view: nglview.NGLWidget
        The view with the 2 molecules
    """
    import nglview
    if radius is None:
        radius = 2.5
    view = nglview.NGLWidget()
    struct1 = view.add_structure(nglview_struct(molecule_low_res))
    struct1.clear_representations()
    struct1.add_representation("licorice", radius=radius, opacity=0.3)

    struct2 = view.add_structure(nglview_struct(molecule_high_res))
    struct2.clear_representations()
    struct2.add_representation("licorice")
    return view

def compare_alignment(manager: Manager, radius: float=None):
    """
    Creates a vertical stack of views for comparing the results of the whole
    alignment.

    Parameters
    ----------
    manager: Manager
        The manager of the alignment.
    radius: Optional[float]
        The radius of the atoms in the low resolution representation. Default
        2.5

    Returns
    -------
    box: ipywidgets.Box
        A vertical box with al the views ready to be visualizaed in jupyter.


    """
    items = []

    for specie, correspondence in manager.complete_correspondence.items():
        if correspondence.end is None or correspondence.start is None:
            continue
        items.append(Label(value=specie))
        items.append(compare_molecules(correspondence.start,
                                         correspondence.end,
                                         radius=radius))
    return VBox(items)
















