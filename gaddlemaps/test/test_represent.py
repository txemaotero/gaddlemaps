import pytest
import os.path

from gaddlemaps import  interactive_restrictions, compare_alignment, Manager
from gaddlemaps.components import Molecule

from ipywidgets import VBox, Tab, Accordion

# Define some file names
ACTUAL_PATH = os.path.split(os.path.join(os.path.abspath(__file__)))[0]
BMIM_AA_ITP = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.itp')
BMIM_AA_GRO = os.path.join(ACTUAL_PATH, '../data/BMIM_AA.gro')

BF4_AA_ITP = os.path.join(ACTUAL_PATH, '../data/BF4_AA.itp')
BF4_AA_GRO = os.path.join(ACTUAL_PATH, '../data/BF4_AA.gro')

BMIM_CG_ITP = os.path.join(ACTUAL_PATH, '../data/BMIM_CG.itp')
BF4_CG_ITP = os.path.join(ACTUAL_PATH, '../data/BF4_CG.itp')

BMIMBF4_SYS = os.path.join(ACTUAL_PATH, '../data/system_bmimbf4_cg.gro')



@pytest.fixture
def manager_bmim():
    manager = Manager.from_files(BMIMBF4_SYS, BMIM_CG_ITP, BF4_CG_ITP)
    bmim = Molecule.from_files(BMIM_AA_GRO, BMIM_AA_ITP)
    bf4 = Molecule.from_files(BF4_AA_GRO, BF4_AA_ITP)

    manager.molecule_correspondence['BMIM'].end = bmim
    manager.molecule_correspondence['BF4'].end = bf4

    return manager

@pytest.fixture
def manager_missing_bmim():
    manager = Manager.from_files(BMIMBF4_SYS, BMIM_CG_ITP, BF4_CG_ITP)
    bf4 = Molecule.from_files(BF4_AA_GRO, BF4_AA_ITP)

    manager.molecule_correspondence['BF4'].end = bf4

    return manager

def test_restrictions_styles(manager_bmim):

    widget, restriction = manager_bmim.interactive_restrictions()

    assert isinstance(widget, VBox)
    assert "BMIM" in restriction
    assert "BF4" in restriction
    assert len(restriction) == 2

    widget, restriction = manager_bmim.interactive_restrictions(style=0)

    assert isinstance(widget, Tab)
    assert "BMIM" in restriction
    assert "BF4" in restriction
    assert len(restriction) == 2

    widget, restriction = manager_bmim.interactive_restrictions(style=1)

    assert isinstance(widget, Accordion)
    assert "BMIM" in restriction
    assert "BF4" in restriction
    assert len(restriction) == 2

    widget, restriction = manager_bmim.interactive_restrictions(style=2)

    assert isinstance(widget, VBox)
    assert "BMIM" in restriction
    assert "BF4" in restriction
    assert len(restriction) == 2

def test_restrictions_missing(manager_missing_bmim):

    widget, restriction = manager_missing_bmim.interactive_restrictions()

    assert isinstance(widget, VBox)
    assert "BMIM" not in restriction
    assert "BF4" in restriction

def test_compare(manager_bmim, manager_missing_bmim):
    compare = compare_alignment(manager_bmim)

    assert len(compare.children) == 4

    compare_missing = compare_alignment(manager_missing_bmim)

    assert len(compare_missing.children) == 2
