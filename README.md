# GADDLE Maps: General Algorithm for Discrete Object Deformations Based on Local Exchange Maps

[**Installation**](#installation)
|
[**Quick Start**](#quick-start)
|
[**Features**](#features)
|
[**FAQs**](#faqs)

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
backend](#configuring-the-c++-backend) first.

The installation can be done using [pip](#pip-installation) or [building from
source](#building-from-source). For any of the
methods gaddlemaps needs [Python](http://pythonhosted.org) 3.6 or greater.

### pip installation

**NOTE: THIS DOES NOT WORK YET, IT HAS NOT BEEN UPLOADED TO PIPY**

This is the simplest method to install the module just running

```bash
pip install gaddlempas --upgrade
```

### Building from source
#### Prerequisites

In order to build from source the following prerequisites are needed:

* [Numpy](https://numpy.org)
* [Scipy](https://www.scipy.org)
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

If it is not installed refer to the section [Configuring the c++ backend](#Configuring-the-c++-backend).

### I have a big protein and it takes ages to align. Is this behavior normal?

First of all check wether you are using or not the c++ backend [How do I know if
gaddlemaps is using the c++
backend?](#How-do-I-know-if-gaddle-maps-is-using-the-c++-backend?). However even
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