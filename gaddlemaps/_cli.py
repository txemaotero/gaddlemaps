import argparse
from typing import Sequence, List, Dict, Set

def auto_map(refrence_coordinates: str, species: Sequence[Sequence[str]],
             scale: float=0.5, outfile: str=None):
    from . import Manager
    from .parsers import read_topology
    from .components import Molecule
    import os.path
    
    folder, basename = os.path.split(refrence_coordinates)
    
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
    
    if outfile is None:
        out_path = os.path.join(folder, f"mapped_{basename}")
    else:
        out_path = outfile
    manager.extrapolate_system(out_path)
    
    return

def classify_files(files):
    from .parsers import ParserManager
    from .parsers._top_parsers import TopologyParserManager
    import os.path
        
    topology_files = set()
    coordinate_files = set()
    
    for filename in files:
        name = os.path.basename(filename)
        extension = name.split(".")[-1]
        
        # Right now there is no file that can be a coordinate and a topology
        # file at the same time, but they may exist in the future or an user may
        # create a custom one, and therefore we allow the same file to be
        # coordinate and topology.
        if extension in ParserManager.parsers:
            coordinate_files.add(filename)
        
        if extension in TopologyParserManager.parsers:
            topology_files.add(filename)
    
    return topology_files, coordinate_files

def sort_molecules(reference_coordinates: str, all_files: Sequence[str],
                   known_files: Sequence[Sequence[str]]) -> Dict[str,Dict[str, str]]:

    from .components import System, MoleculeTop, Molecule
    import warnings
    
    # Let's clasify the files in topology and coordinates
    topology_files, coordinate_files = classify_files(all_files)

    
    # Remove the known files to avoid duplicates
    for molecule_files in known_files:
        if molecule_files[0] in topology_files:
            topology_files.remove(molecule_files[0])

        if molecule_files[1] in coordinate_files:
            coordinate_files.remove(molecule_files[1])
        
        if molecule_files[2] in topology_files:
            topology_files.remove(molecule_files[2])
    
    # Initiate the system and add the known topologies
    
    system = System(reference_coordinates, *[files[0] for files in known_files])
    
    added_molecues: Dict[str,Dict[str, str]] = {}
    used_files: Set[str]  = set()
    
    topology_molecues = [(i, MoleculeTop(i)) for i in topology_files]
    
    for filename, molecule in topology_molecues:
        try:
            system.add_molecule_top(molecule)
        except OSError:
            pass
        else:
            used_files.add(filename)
            added_molecues[molecule.name] = {"top_CG": filename}
            
    for filename, molecule in topology_molecues:
        if (filename not in used_files) and (molecule.name in added_molecues):
            if "top_AA" in added_molecues[molecule.name]:
                warnings.warn(f"Repeated topology found for molecule {molecule.name}",
                              RuntimeWarning)
            else:
                added_molecues[molecule.name]["top_AA"] = filename
    
    # Look for a matching coordinate file
    
    for coordinate_file in coordinate_files:
        for molecule_name, molecule_info in added_molecues.items():
            if "coor_AA" not in molecule_info:
                try:
                    Molecule.from_files(coordinate_file, molecule_info["top_AA"])
                except OSError:
                    pass
                else:
                    added_molecues[molecule_name]["coor_AA"] = coordinate_file
    
    return added_molecues
    
    
    
    
    
    
    
    
    

def main():
    description = """
    Fill me
    """
    
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("init_coor", type=str, metavar="init_coordinates.gro",
                        help=("Add doc here"))
    
    parser.add_argument('-mol', dest="mol", action='append', nargs=3, type=str,
                        metavar="mol_cg.itp mol_aa.gro mol_aa.itp")
    
    parser.add_argument('--scale', dest="scale", default=0.5, type=float,
                        metavar="0.5")
    
    parser.add_argument("-o", "--outfile", dest="outfile", required=False)
    
    parser.add_argument("--auto", dest="auto", action="store", nargs="+", type=str)
    
    
    args = parser.parse_args()
    
    if args.mol is None:
        molecules = []
    else:
        molecules = args.mol
    
    if args.auto is not None:
        print("Starting automatic search of molecues:\n")
        molecule_info = sort_molecules(args.init_coor, args.auto, molecules)
        
        count_mols = sum([len(i) == 3 for i in molecule_info.values()])
        print(f"A total of {count_mols} molecules has been automatically added\n")
        
        for molecule_name in molecule_info:
            if len(molecule_info[molecule_name]) == 3:
                new_mol = [molecule_info[molecule_name]["top_CG"],
                           molecule_info[molecule_name]["coor_AA"],
                           molecule_info[molecule_name]["top_AA"]]
                molecules.append(new_mol)
                text = ("The molecue {} has been added with initial topology {}"
                        ", final coordinates {} and final topology {}")
                print(text.format(molecule_name, *new_mol))
        print("\n\n")
                
    auto_map(args.init_coor, molecules, args.scale, outfile=args.outfile)