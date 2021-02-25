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
import argparse
from typing import Sequence, Dict, Set, Tuple

def auto_map(refrence_coordinates: str, species: Sequence[Sequence[str]],
             scale: float=0.5, outfile: str=None):
    """
    Performs all the steps in the mapping process
    """
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

def classify_files(files) -> Tuple[Set[str], Set[str]]:
    """
    Clasify the files betweeen topology and coordinate files

    Paramters
    ---------
    files: Sequence[str]
        All the files that will be clasified

    Returns
    -------
    topology_files: Set[str]
        The files compatible with a topology parser

    coordinate_files: Set[str]
        The files compatible with a coordinate parser
    """
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
    """
    Given a list of unknown files looks for molecules compatible with the
    reference coordiantes and that are not included in the knwon files
    """

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

    added_molecues: Dict[str, Dict[str, str]] = {}
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
    from gaddlemaps import check_backend_installed

    description = """
    Command line interface for the python module gaddlemaps.

GADDLE Maps(https://doi.org/10.1021/acs.jctc.7b00861) is a general algorithm
 to switch between structures with different resolutions (*i.e.* different number
of particles per object). This algorithm is completely automatic and does not
require any human input. The algorithm has been designed specifically to
transform molecules between to levels of resolution (for example from united
atom to all atom or from coarse grain to fully atomistic). However the algorithm
works for any structure than can be represented by a graph.

This implementation allows an easy mapping between molecular systems in two
different resolutions. The command line interface is the easiest way to map the
simplest systems. For more options and the ability to graphically select the
molecule restraints please refer to the README and the examples provided with
the source code.

The command line interface has also the ability to automatically select the
topology and coordinate files necessary to perform the mapping from a list of
possible files.

    """

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("init_coor", type=str, metavar="init_coordinates.gro",
                        help=("The inital system in the starting resolution."
                        "This filee must contain all the coordinates of the "
                        "system"))

    help_mol = """Files with the information of each
    molecule. Must provide 3 files in the following order: topology file in the
    starting resolution, coordinate file in the final resolution and the
    topology file in the final resolution. This argument must be called once for
    each specie in the system that you want to map. For example '--mol
    BMIM_CG.itp BMIM_AA.gro BMIM_AA.itp --mol BF4_CG.itp BF4_AA.gro BF4_AA.itp'"""

    parser.add_argument('-m', "--mol", dest="mol", action='append', nargs=3, type=str,
                        metavar=("mol_cg.itp", "mol_aa.gro", "mol_aa.itp"),
                        help=help_mol)

    help_scale = """Scaling factor that shrinks the system (for values <1).
    This shrink factor is used to avoid colision between the atoms in the final
    resolution. It is usually needed because the volume that the molecule fills
    in both resolutions is not usually the same. A short energy minimization in
    the final resolution is usually enough to recover a "normal" system."""

    parser.add_argument('--scale', dest="scale", default=0.5, type=float,
                        metavar="0.5", help=help_scale)

    help_outfile = """Output file where the mapped system will be saved.
    Default mapped_$initial_file"""

    parser.add_argument("-o", "--outfile", dest="outfile", required=False,
                        help=help_outfile)

    help_auto = """A list of files. gaddlemaps will scan this list to find
    information about the molecules in the system. This works only if the name
    of the molecule in topology files is the same in both resolutions. This
    option can be combined with the -mol option. Files for the molecules set
    with the -mol flag will not be searched."""

    parser.add_argument("--auto", dest="auto", action="store", nargs="+",
                        type=str, help=help_auto)

    help_exclude = """Molecule names of species that you don't want to map. This
     flag works with the --auto flag. If gaddlemaps finds a molecule in the
     files given in the --auto flag but that molecule is included in the exclude
      flag it will not be added to the mapping. This is usually useful to avoid
      mapping solvent molecules (for example watter) that can be easily readded
      after the mapping."""

    parser.add_argument("--exclude", dest="exclude", action="store", nargs="+",
                        type=str, help=help_exclude)


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

                if (args.exclude is not None) and (molecule_name in args.exclude):
                    print(f"Excluding automatically found molecule {molecule_name}")
                    continue

                new_mol = [molecule_info[molecule_name]["top_CG"],
                           molecule_info[molecule_name]["coor_AA"],
                           molecule_info[molecule_name]["top_AA"]]
                molecules.append(new_mol)
                text = ("The molecue {} has been added with initial topology {}"
                        ", final coordinates {} and final topology {}")
                print(text.format(molecule_name, *new_mol))
        print("\n\n")

    backend = check_backend_installed(warn_missing=True)

    if backend:
        print("Starting alingment using optimized backend")
    else:
        print("Optimized backend was not found, falling back to the python implementation.")


    auto_map(args.init_coor, molecules, args.scale, outfile=args.outfile)