# GADDLE Maps: General Algorithm for Discrete Object Deformations Based on Local Exchange Maps

[**Installation**](#installation)
|
[**Quick Start**](#quick-start)
|
[**Features**](#features)
|
[**FAQs**](#faqs)
|
[**Documentation**](https://gaddlemaps.readthedocs.io/)

## What is GADDLE Maps?

[GADDLE Maps](https://doi.org/10.1021/acs.jctc.7b00861) is a general algorithm
to switch between structures with different resolutions (*i.e.* different number
of particles per object). This algorithm is completely automatic and does not
require any human input. The algorithm has been designed specifically to
transform molecules between to levels of resolution (for example from united
atom to all atom or from coarse grain to fully atomistic). However the algorithm
works for any structure than can be represented by a graph.

The algorithm apply rotations and deformations to the graph (molecule) preserving the
distance between its nodes (atoms) in order to maximize the overlapping between
the graphs (molecules) in both resolutions. Once this is achieved, it generates
an *exchange map* that transforms between both resolutions. These exchange maps
use local information about the location of the neighbours to exchange between
the 2 resolutions. This allows the algorithm to adapt to multiple instances of
the same graph (molecule) in different configurations. Using this property the
algorithm is able to transform subsequent instances of the same graph (molecule)
at almost zero cost. For these reason the algorithm is specially suitable to
transform between systems in coarse grained and fully atomistic resolutions
because by only aligning one molecule of each specie the whole system can be
transformed.

> **IMPORTANT** If you are experiencing problems during the mapping process (the
> mapping gets stuck at some point, the final system looks weird...) check the
> [FAQs](#FAQs) section because most common mistakes are covered there.

### The gaddlemaps module

In this module we offer the reference implementation of the algorithm in python
with a series of helper classes and functions for an easy usage of the
algorithm. Moreover a command line interface and a jupyter notebook is provided
for an easy and visual way of performing the transformation.

Besides a port of the algorithm to c++ with an interface with python using
[cython](https://cython.org). This port allows for a boosted performance,
specially for large molecules such as proteins, but the results are the same
than as the reference implementation. This port can be installed alongside
the module. `gaddlemaps` will automatically use the c++ port if available and will
fallback to the python implementation if not present. Therefore the installation
of this c++ port is completely **optional** but highly recommended.

### Content
* [Installation](#installation)
* [Quick Start](#quick-start)
* [Frequently asked questions](#FAQs)

## Installation

Here we will describe the installation of the reference implementation of
`gaddlemaps` (*i.e.* without the c++ backend). In order to install the
backend please read the section [configuring the c++
backend](#configuring-the-c-backend) first.

The installation can be done using [pip](#pip-installation) or [building from
source](#building-from-source). For any of the
methods gaddlemaps needs [Python](http://pythonhosted.org) 3.6 or greater.

### pip installation

This is the simplest and prefered method to install the module. In order to
install gaddlemaps just run:

```bash
pip install gaddlempas --upgrade
```

This method will autodetect wether the requirements for the optimized c++
backend are available and, if possible, will install the module with the
optimized backend. If the dependencies are not met, it will install the native
python version of the backend. The dependencies necessary for the backend are
detailed in the section [configuring the c++
backend](#configuring-the-c-backend).

### Building from source
#### Prerequisites

In order to build from source the following prerequisites are needed:

* [Numpy](https://numpy.org)
* [Scipy](https://www.scipy.org)
* [more_itertools](https://github.com/more-itertools/more-itertools)
* [Jupyter Notebook](https://jupyter.org)
* [nglview](https://github.com/arose/nglview) (optional)

#### Installation of prerequisites

This requisites can installed using **pip** by running:

```bash
pip install numpy scipy jupyter
```

If using **Mac OS** the dependencies can be installed using [brew](https://brew.sh) by
running:

```bash
brew install numpy scipy jupyter
```

In **Ubuntu** and other Debian based Linux distribution they can be installed using
the *apt* package manager:

```bash
sudo apt install python3-numpy python3-scipy python3-notebook
```

If using a different operative system or if the previous instructions did not
work please refer to the installation instructions in the documentation of each
module.

**Nglview** is completely optional, it is used only when representing the system
inside jupyter notebooks. If you do not want this functionality you can skip
this step. Moreover this module **can be installed after installing
gaddlemaps**. The easiest way to install the module is through pip

```bash
pip install nglview
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter-nbextension enable nglview --py --sys-prefix
```

for other options and problems with the installation check the documentation for
[nglview](https://github.com/arose/nglview).

#### Cloning and installing

The first step is cloning or downloading the git repository

```bash
git clone https://github.com/txemaotero/gaddlemaps
cd gaddlemaps
```

or

```bash
wget https://github.com/txemaotero/gaddlemaps/archive/master.zip
unzip master.zip
cd gaddlemaps-master
```

After that the module is installed by running:

```bash
sudo python3 setup.py install
```

The module can also be installed without admin privileges with:

```bash
python3 setup.py install --user
```

After that the module will be installed and the downloaded folder can be safely
removed.

### Configuring the c++ backend

In order to install the c++ backend the following extra prerequisites are
needed:

* A c++ compiler
* [Armadillo](http://arma.sourceforge.net)
* [Cython](https://cython.org)

#### C++ compiler

The installation of a c++ compiler (for example g++ or clang++) is needed. This
can be done in **Ubuntu** and other Debian-like Linux distribution using the apt package
manager:

```bash
sudo apt install g++
```

If using **Mac OS** the compiler can be installed by installing the *developer command
line tools*. These tools can be installed by running the following command in a
terminal:

```bash
xcode-select --install
```

and following the instructions on the screen.

#### Armadillo

For the installation of the Armadillo library please refer to their
[installation documentation](http://arma.sourceforge.net/download.html).

If using a Debian based Linux distribution it should be possible to install it
using apt:

```bash
sudo apt-get install libarmadillo-dev
```

In order to check if the installation is in the path create a file
*armadillo_test.cpp* with the following content:

```c++
#include <armadillo>
#include <iostream>

int main(){
    arma::Mat<double> A(3,3, arma::fill::eye);
    std::cout << "A:\n" << A << "\n";
    return 0;
}
```

And try to compile and run it:

```bash
g++ armadillo_test.cpp -larmadillo -I/usr/local/include -L/usr/local/lib && ./a.out
```

If a 3⨉3 identity matrix is printed the armadillo library is installed and in a
location where gaddlemaps will be able to find it. If the compiler is not able
to find the armadillo library update the environment variables in your
*.bashrc* (if using Linux or Mac OS prior to Catalina) or *.zshrc* (if using zsh instead of bash or if using Mac OS Catalina).

In general adding the following lines at the end of the file will de the trick:

```bash
export CPLUS_INCLUDE_PATH="folder where armadillo.h is located"
export LD_LIBRARY_PATH="folder where libarmadillo.so  or libarmadillo.dylib"
export LIBRARY_PATH="folder where libarmadillo.so  or libarmadillo.dylib"
```

If armadillo was installed in Mac OS using [macports](https://www.macports.org)
the following lines should work:

```bash
export CPLUS_INCLUDE_PATH="/opt/local/include"
export LIBRARY_PATH="/opt/local/lib"
export LD_LIBRARY_PATH="/opt/local/lib"
```

#### Cython

##### Using pip

The easiest and preferred way to install cython is by using pip

```bash
sudo pip install cython
```

##### Without pip

In **Mac OS** can be installed using brew

```bash
brew install cython
```

In **Ubuntu** can be installed using apt

```bash
sudo apt install cython3
```

For other platforms and other installation method check the official
documentation for [cython](https://cython.org).

After finishing the installation of the prerequisites please follow the steps in
the [installation](#installation) section.

## Quick Start

### Files needed to map the system

gaddlemaps uses 2 different kind of files *coordinates files* and *topology
files*. The coordinates file stores the coordinates of the atoms while the
topology file stores the connections (*i.e.* bonds) between the atoms. Currently
the topology file can only store the information for 1 molecule, while the
coordinate file can store the whole system.

For doing the whole mapping process the following information is required:

- The coordinate file with the whole system in the start resolution.
- One topology file for each specie that you want to map (it is not required to
  map all the molecules in the system).
- For each specie that you want to map a coordinate file with the molecule in
  the end resolution (it is possible to use only 1 coordinate file with all the
  molecules, but it requires more intervention by the user).
- One topology file for each specie that you want to map in the end resolution.

gaddlemaps has been developed using [GROMACS](http://www.gromacs.org) as
reference molecular dynamics simulation suite. Therefore currently the allowed
formats for the coordinate and topology files are:

- Coordinate Files: *.gro*
- Topology files: *itp*

However, it is easy to add more formats without touching any of the code in the
module by the user. If you are interested in adding support for other formats
read the FAQ [How can I use an unsupported file format?](#how-can-i-use-an-unsupported-file-format)

### Perform the mapping

The basic usage of the module to transform systems between coarse grained and
fully atomistic resolutions can be seen in the example
[examples/bmim_bf4/map_bmim_bf4.py](examples/bmim_bf4/map_bmim_bf4.py).

For a more visual representation there is a jupyter notebook that works as
example and can be easily modified to transform interactively any system. The
notebook is located in
[examples/notebook/Interactive-example.ipynb](examples/notebook/Interactive-example.ipynb).
In order to use this notebook the [nglview](https://github.com/arose/nglview)
module is required. In order to run the notebook run the following command:

```bash
jupyter notebook examples/notebook/Interactive-example.ipynb --NotebookApp.iopub_data_rate_limit=10000000
```

If you do not have [Jupyter](https://jupyter.org) installed you can take a look
to its markdown preview in
[examples/notebook/Interactive-example.md](examples/notebook/Interactive-example.md).

If you are looking for a more flexible way of performing the mappings that
allows you to have more control in each step you can find another Jupyter notebook
[example](examples/notebook/manual_mapping.ipynb) (and the corresponding
[preview](examples/notebook/manual_mapping.md)) which will take you around the
module to explore more options during the mapping process. This example will
also tell you how you can reuse the results from old mappings to afford
computational time.

#### CLI

Another option to perform a mapping is using the command line interface that is
available after the module installation. This tool is accessible through the
command `gaddlemaps` in a terminal. The documentation of this command can be
checked running:

```bash
gaddlemaps -h
```

A basic usage of it would look like this:

```bash
gaddlemaps system_initial.gro --mol molecule1_initial.itp molecule1_final.gro molecule1_final.itp --mol molecule2_initial.itp molecule2_final.gro molecule2_final.itp --scale 0.5 -o system_final.gro
```

where `system_initial.gro` is the file with the coordinates of the atoms of the
system to be mapped in the initial resolution. After each `--mol` flag all the
files for mapping each desired type of molecule must be provided (topology in
the initial resolution, atom coordinates in the final resolution and topology of
the molecule in the final resolution). The `--scale` flag specifies the scaling
factor to reduce the final atomic coordinates to avoid molecular overlaps.
Finally, the `-o` option specified the path to write the output file with the
coordinates in the final resolution.

There is another way to run a mapping with this tool that automatically detects
all the files needed for each molecule. For example, imagine that in the current
directory you have all the files needed for the mapping we did in with the
previous command and more (for instances files to map a "molecule3" that we do
not want to map). In this scenario we could run the same mapping as before
running the following command:

```bash
gaddlemaps system_initial.gro --auto * --exclude molecule3 --scale 0.5 -o system_final.gro
```

Note that the identification of the molecules is based on their names in the
".itp" files and also in the possibility of molecule object initiation so use
this option with caution. The `--exclude` flag is used to tell the command to
not include the molecule3 (the molecule name in the .itp file) in the mapping
process (useful for example to exclude solvent molecules).


## FAQs

### What are those restraints? Wasn't gaddlemaps an automatic algorithm?

GADDLE Maps is an algorithm that does not need any external input. However for
molecules that are very symmetrical (for example molecules that are completely
linear) the algorithm may need a "hint" to know which extreme is which. When
the molecule have ramification these hints are not usually needed. These hints
are what we call restraints. They can be also used in order to have a faster aligment.

### How do I know if gaddle maps is using the c++ backend?

There is a function that returns wether the c++ backend is available or not. In
order to check it run:

```python
from gaddlemaps import check_backend_installed
if check_backend_installed():
    print("The backend is correctly installed")
else:
    print("The backend is not accessible")
```

If it is not installed refer to the section [Configuring the c++ backend](#Configuring-the-c-backend).

### I have a big protein and it takes ages to align. Is this behavior normal?

First of all check wether you are using or not the c++ backend [How do I know if
gaddlemaps is using the c++
backend?](#How-do-I-know-if-gaddle-maps-is-using-the-c-backend). However even
without the c++ backend the mapping should not take more than a day even for
really big molecules.

The time that takes to align a molecule greatly increase with the number of
atoms. However, for proteins an approximation is done. Instead of aligning all
the atoms of the protein the overlapping is only calculated between atoms in
residues  with the same name. If the alignment is taken too long chances are
that the residues in both resolutions have different names or are in different
order.

This problem can be solved by adding restraints between all the atoms of the
molecule but by far the easiest way to accelerate the alignment is to set the
same sequence of residues in both files.

### Why the final mapped system looks so weird?

By default, when the molecular alignments are extrapolated to all the
configurations in the simulation box, the coordinates of the atoms in the
molecules in the final resolution are scaled down (see the original paper for
a more detailed description of this reduction) to avoid different molecules
to overlap in the final system. If this is not done, every simulation that
takes that mapped system as input would crash in the first steps. A short energy
minimization is enough to fix this scale of the atomic coordinates.

> **NOTE**: Before trying to resolvate or ionize any of the systems obtained
> with this tool is highly recommended to perform a short energy minimization to
> relax the molecular structures.

### I have some other problem or an unexpected behavior of the code. How should I report it?

If you find that the code does not work as expected please check if there is an
[issue in github](https://github.com/txemaotero/gaddlemaps/issues) with the same
problem. If there is not, [open a new
issue](https://github.com/txemaotero/gaddlemaps/issues/new/choose) and we will
address it as soon as possible.

### How can I use an unsupported file format?

gaddlemaps allows the enduser to extend the compatible file formats without
touching a single line of the module. For such purpose gaddlemaps have 2 base
classes (one for the coordinate files and another for the topology files) from
which the enduser can inherit. Any class that inherit from these classes will be
added automatically to the supported formats.

Any class that inherits from *gaddlemaps.parsers.CoordinatesParser* will be
added to the compatible parsers for coordinate files. On the other any class
that inherits from *gaddlemaps.parsers.TopolgyParser* will be added to the
supported parsers for topology files. Below you can see an example for a simple
coordinate parser for the *pdb* format (the *pdb* format is not currently
supported as the parser below is not guarantee to work for the whole [pdb
specification](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)).

**NOTE:** The following snippet is just a quick and dirty implementation of a
parser. It is just an example of how a quick parser can be written. In its
current form quality would not be good enough to be included as part of the
module.

```python
from gaddlemaps.parsers import CoordinatesParser

import numpy as np

class SimplePDB(CoordinatesParser):

    # We must define wich extensions will the parser work with
    EXTENSIONS = ("pdb", )

    # The following methods must be defined

    def __init__(self, path, mode='r'):
        self.open_file = open(path, mode=mode)

        self.mode = mode
        self._comment = ""
        self._box = np.array([[0, 0, 0],
                              [0, 0, 0],
                              [0, 0, 0]])
        self._atom_info = []
        self._atom_index = 0

        if mode == "r":
            self._load_info()

    def seek_atom(self, index):
        # Instead of moving the cursor to the desired atom we can read all the
        # atomic information first and use an index to seek the atom. Beware
        # that this will increase the needed ram and can give problems for the
        # bigger systems.

        self._atom_index = index

    def next(self):
        if self._atom_index >= self.natoms:
            raise StopIteration

        info = self._atom_info[self._atom_index]

        self._atom_index += 1
        return info

    def writeline(self, atomlist):
        # Instead of actually writing the file on the go it is possible to
        # store the information and save it all just before closing

        self._atom_info.append(atomlist)



    def close(self):
        if self.mode == "w":
            self.open_file.write(f"TITLE   {self.comment}\n")
            self.open_file.write("REMARK    THIS IS A SIMULATION BOX\n")
            self.open_file.write("MODEL 1\n")
            for atomlist in self._atom_info:
                atom_index = atomlist[3]
                atom_name = atomlist[2]
                residue_name = atomlist[1]
                residue_index = atomlist[0]
                x_position = atomlist[4] * 10  # The original value is in nm
                y_position= atomlist[5] * 10
                z_position = atomlist[6] * 10


                line = (f"ATOM {atom_index:6d} {atom_name} {residue_name} {residue_index:6d}"
                        " {x_position:8.3f} {y_position:8.3f} {z_position:8.3f}\n")
                self.open_file.write(line)
            self.open_file.write("TER\n")
            self.open_file.write("ENDMDL\n")

        self.open_file.close()

    # We must define some properties also

    @property
    def natoms(self):
        return len(self._atom_info)

    @property
    def box_matrix(self):
        return self._box

    @box_matrix.setter
    def box_matrix(self, new_box):
        self._box = np.array(new_box)

    @property
    def comment(self):
        return self._comment

    @comment.setter
    def comment(self, newcomment):
        self._comment = newcomment

    # You can add as many auxiliar methods as you want

    def _load_info(self):
        for line in self.open_file:
            line = line.strip()
            if line.startswith("TITLE"):
                self.comment = line[5:].strip()

            elif line.startswith("ATOM"):
                line = line.split()
                self._atom_info.append((
                    int(line[4]),
                    line[3],
                    line[2],
                    int(line[1]),
                    float(line[5])/10,
                    float(line[6])/10,
                    float(line[7])/10
                    ))
            else:
                pass
```

Just having this class imported or *pasted* in your script will allow gaddlemaps
to be used with *.pdb* coordinate files.