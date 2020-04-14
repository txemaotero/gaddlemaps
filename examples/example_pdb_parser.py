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
                    line[1],
                    float(line[5])/10,
                    float(line[6])/10,
                    float(line[7])/10
                    ))
            else:
                pass
    
test = SimplePDB("/usr/local/il-library/opls_aa/cation-aprotic/BMIM.pdb")

for atom in test:
    print(atom)