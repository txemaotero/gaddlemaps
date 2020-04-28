import argparse
from typing import Sequence

def auto_map(refrence_coordinates: str, species: Sequence[Sequence[str]],
             scale: float=0.5):
    from . import Manager
    from .parsers import read_topology
    from .components import Molecule
    
    # Load the molecule names
    end_molecules = {}
    itps_cg = []
    
    for specie in species:
        name, *_ = read_topology(specie[0])
        
        end_molecules[name] = Molecule.from_files(specie[1], specie[2])
        itps_cg.append(specie[0])
    
    manager = Manager.from_files(refrence_coordinates, *itps_cg)
    
    for name, end_mol in end_molecules.items():
        manager.molecule_correspondence[name].end = end_mol
    
    manager.align_molecules()
    
    manager.calculate_exchange_maps(scale_factor=scale)
    
    manager.extrapolate_system(f"{refrence_coordinates}_mapped.gro")
    
    return

def main():
    description = """
    Fill me
    """
    
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("init_coor", type=str, metavar="init_coordinates.gro",
                        help=("Add doc here"))
    
    parser.add_argument('-mol', dest="mol", action='append', nargs='+', type=str,
                        metavar="mol_cg.itp mol_aa.gro mol_aa.itp")
    
    parser.add_argument('--scale', dest="scale", default=0.5, type=float,
                        metavar="0.5")
    
    
    args = parser.parse_args()
    
    auto_map(args.init_coor, args.mol, args.scale)