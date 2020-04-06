# GADDLE Maps: General Algorithm for Discrete Object Deformations Based on Local Exchange Maps

[**Installation**](#installation)
|
[**Quick Start**](#quick-start)
|
[**Features**](#features)

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

Moreover a port of the algorithm to c++ with an interface with python using
[cython](https://cython.org). This port allows for an boosted performance,
specially for large molecules such as proteins, but the results are the same
than using the reference implementation. This  port can be installed alongside
the module. gaddlemaps will automatically use the c++ port if available and will
fallback to the python implementation if not resent. Therefore the installation
fo this c++ port is completely **optional**.

## Installation

Here we will describe the installation of the reference implementation of
gaddlemaps (*i.e.* without the c++ backend). In order to install the backend
please read the section [configuring the c++
backend](#configuring-the-c++-backend) first.

The installation can be done using [pip](#pip-installation) or building from source. For any of the
methods gaddlemaps needs [Python](http://pythonhosted.org) 3.6 or greater.

### pip installation

**NOTE: THIS DOES NOT WORK YET, IT HAS NOT BEEN UPLOADED TO PIPY**

This is the simplest method to install the module just running

```bash
pip install gaddlempas --upgrade
```

### Building from source
#### Prerequisites

In order to build from source the following prerequisites are needeed:

* [Numpy](https://numpy.org)
* [Scipy](https://www.scipy.org)
* [Jupyter Notebook](https://jupyter.org)

This requisites can installed using pip by running:

```bash
pip install numpy scipy jupyter
```

If using mac the dependencies can be installed using [brew](https://brew.sh) by
running:

```bash
brew install numpy scipy jupyter
```

In Ubuntu and other Debian based Linux distribution they can be installed using
the *apt* package manager:

```bash
sudo apt install python3-numpy python3-scipy python3-notebook
```

If using a different operative system or if the previous instructions did not
work please refer to the installation instructions in the documentation of each
module.

### Cloning and installing

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

The installation of a c++ compiler (for example g++ or clang++) is needed. This
can be done in Ubuntu and other Debian-like Linux distribution using the apt package
manager:

```bash
sudo apt install g++
```

If using Mac the compiler can be installed by installing the *developer command
line tools*. These tools can be installed by running the following command in a
terminal:

```bash
xcode-select --install
```

and following the on screen instructions.

For the installation of the Armadillo library please refer to their
[installation documentation](http://arma.sourceforge.net/download.html)

In order to check if the installation is in the path create a file
*armadillo_test.cpp* with the folllowing content:

```c++
#include <armadillo>
#include <iostream>

int main(){
    arma::Mat<double> A(4,4, arma::fill::eye);
    std::cout << "A:\n" << A << "\n";
    return 0;
}
```

And try to compile and run it:

```bash
g++ armadillo_test.cpp -larmadillo && ./a.out
```