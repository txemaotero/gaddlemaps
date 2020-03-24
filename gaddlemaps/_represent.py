from . import Manager
from ._manager import Restrictions
from typing import List, Tuple, Dict
from ipywidgets import Button, HBox, VBox, Layout, Label, Box, Tab, Accordion, Widget
from os.path import join, isdir
from os import mkdir

Restriction = List[Tuple[int, int]]



class Cycle:
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

def create_widget_restrictions(file_low_res: str, file_high_res: str,
                                restrict: Restriction = None)-> Tuple[Box, Restriction]:
    import nglview

    if restrict is not None:
        restrictions = restrict
    else:
        restrictions = []
    
    
    
    text = Label(value="Low Resolution")
    view = nglview.show_file(file_low_res)

    text2 = Label(value="High Resolution")
    view2 = nglview.show_file(file_high_res)

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
                color=next(colors)
                view.add_representation("licorice", radius=1, color=color, selection=[restriction[0]], opacity=0.75)
                view2.add_representation("licorice", radius=1, color=color, selection=[restriction[1]], opacity=0.75)


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


def create_interactive_restriction(manager: Manager,
                                   folder: str=None) -> Tuple[Dict[str, Box],
                                                              Restrictions]:

    if folder is None:
        folder = "representation_tempfiles"
    
    if not isdir(folder):
        mkdir(folder)
    
    restrictions = {}
    boxes = {}
    
    for specie in manager.molecule_correspondence:
        if manager.molecule_correspondence[specie].end is None:
            continue
        low_res_file = join(folder, f"{specie}_low_res.gro")
        high_res_file = join(folder, f"{specie}_high_res.gro")
        
        manager.molecule_correspondence[specie].start.write_gro(low_res_file)
        manager.molecule_correspondence[specie].end.write_gro(high_res_file)
        
        box, restrict = create_widget_restrictions(low_res_file, high_res_file)
        
        boxes[specie] = box
        restrictions[specie] = restrict
        
    
    
    return boxes, restrictions

def interactive_restrictions(manager: Manager, folder: str=None,
                             style:int=2) -> Tuple[Widget, Restrictions]:
    
    if style == 0:
        representation = Tab
    elif style == 2:
        representation = VBox
    else:
        representation = Accordion
    
    boxes, restrictions = create_interactive_restriction(manager, folder)
    
    restriction_widget = representation(sync=False)
    
    for index, specie in enumerate(boxes):
        child = boxes[specie]
        if style in [2, ]:
            child.children = (Label(value=r"$\textbf{"+ specie + r"}$"), ) + child.children
        restriction_widget.children += (boxes[specie],)
        if style not in [2, ]:
            restriction_widget.set_title(index, specie)
    
    return restriction_widget, restrictions
    
    
    
    
    
        
        
        
    
    



