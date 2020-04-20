In this example, a simple simulation system with 300 ionic pairs of [BMIM][BF4]
ionic liquid is mapped from coarse-grained (CG) to fully atomistic (AA)
resolution. Here, the Manager class is used to achieve an easier implementation
of the mapping.


```python
from gaddlemaps import Manager, interactive_restrictions, compare_alignment
from gaddlemaps.components import Molecule
```

First, a Manager instance needs to be created. The easiest way to do that is
using the "from_files" class method. The first argument of this method is the
path to the gro file with the CG configuration of the system that we want to
map. The following arguments are all the files with the topology information
of the molecules that you want to map. For example, if you have a system with
a lot of water molecules, our advise is not to map these water molecules and
resolvate the system once you have mapped the rest of the molecules. In that
case, you should not pass the topology file of the water to the Manager
initialization.


```python
manager = Manager.from_files(
    '../../gaddlemaps/data/system_bmimbf4_cg.gro',
    '../../gaddlemaps/data/BMIM_CG.itp',
    '../../gaddlemaps/data/BF4_CG.itp'
)
```

Now, we need to specify the molecules in the final resolution (AA). In this
example, we initialize two instances of Molecule, one for the BMIM molecule
and other for the BF4 one. For this task, we create a System object with the
gro and the corresponding topology files. Then, we take the first element of
these systems which are the desired molecules (as the input gro files only
contain one molecule).
We could also have created these molecules with the Molecule.from_file method


```python
bmim_aa = Molecule.from_files('../../gaddlemaps/data/BMIM_AA.gro',
                              '../../gaddlemaps/data/BMIM_AA.itp')
bf4_aa = Molecule.from_files('../../gaddlemaps/data/BF4_AA.gro',
                             '../../gaddlemaps/data/BF4_AA.itp')
```

Once the molecules are created, we need to specify which one correspond with
the molecules in the initial resolution. For that, we have 2 options. First,
we can use the add_end_molecules method of Manager, which will automatically
detect the corresponding molecules based on their names attribute.


```python
manager.add_end_molecules(bmim_aa, bf4_aa)
```

The second option is to assign the molecules manually. For this, we have to
access the element with the correct name in molecule_correspondence attribute
and set the "end" attribute. This option has the advantage that it is not
necessary that the molecules have the same name in both resolution.


```python
manager.molecule_correspondence['BMIM'].end = bmim_aa
manager.molecule_correspondence['BF4'].end = bf4_aa
```

We can add some restrictions to help mapping molecules that have some simmetry.
This restrictions will be used as hints for alignning the molecules.
The easiest way to create the restrictions is using the interactive widget. The usage is very simple,
just click in one of the atoms of the molecule CG and one of the AA molecule and press
the button "add_restriction". When finished move to the next cell.


```python
restriction_widget, restrictions = interactive_restrictions(manager)
restriction_widget
```


    _ColormakerRegistry()



    VBox(children=(VBox(children=(Label(value='$\\textbf{BMIM}$'), Label(value='Low Resolution'), NGLWidget(), Lab…


Now we can perform the aligmnet with the restrictions that we have just created.


```python
manager.align_molecules(restrictions=restrictions)
```

    Aligning BMIM:
    
    Aligning BF4:
    


The execution of the alignment usually is the most time spending step
(although it is not for this example as the molecules are very small).
Once the molecules are aligned, the exchange maps need to be initialized.

The "scale_factor" parameter is used to reduce the final size of the mapped
molecules (if the factor is smaller than 1). This avoids possible molecule
overlapping. However, take into account that this will produce molecular
deformations. This deformations are easily removed with a short energy
minimization.
Finally, we can extrapolate the found alignment to the rest of the molecules
in the system calling the "extrapolate_system" method with the desired output
file name as parameter.


```python
manager.calculate_exchange_maps(scale_factor=0.5)
manager.extrapolate_system('example_bmim_bf4_map.gro')
```

It is very recommended to check the found alignment for each molecule type.
For this task, we can call the "write_comparative_gro" method of the Alignment
object for each molecule. By default, this method will write a file for each
molecule with the name "{molecule.name}_compare.gro". In these gro files you
will find two molecules with the residue names and ids "START", 1 and "END", 2
for the molecules in the initial and final resolution respectively. You can
use any external software of molecular visualization to check if the "END"
molecule correctly overlaps the "START" one. 


```python
for align in manager.molecule_correspondence.values():
    align.write_comparative_gro()
```


```python
compare_alignment(manager, radius=2.5)
```


    VBox(children=(Label(value='BMIM'), NGLWidget(), Label(value='BF4'), NGLWidget()))

