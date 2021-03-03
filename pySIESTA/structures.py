"""
    Provides data structures to handle SIESTA files.
"""

# third party imports
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODULE STRUCTURE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
# + class Geometry()
#   - reset()
#   - read_STRUCT(struct_file)
#   - write_XYZ(xyz_file)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class Geometry():
   
    """
    SIESTA Geometry Container

    Attributes:
    ----------

     - atoms (dict): dictionary with positions of the atoms (fractional)
     - lat_vecs (array): supercell shape 

    """

    def reset():

        self.lat_vecs  = np.zeros((3,3))
        self.atoms     = []
        self.loaded    = False

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def read_STRUCT(struct_file, mode="molecule"):

        self.reset()

        self.mode = mode

        with open(fname, "r") as f:
            
            for i in range(3):
                line = f.readline().split()
                lat_vec[i,:] = [float(x) for x in line]

            self.nats = int(f.readline().split()[0])

            self.pos = np.zeros((2,2,2,5,3))
            
            for at in range(nats):
                            
                line = f.readline().split()

                label   = line[0]
                species = line[1]
                pos = [float(x) for x in line[2:]]

                self.atoms[at] = {
                                    "pos": pos,
                                    "label": label,
                                    "species": species
                                 }
        self.loaded = True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def write_XYZ(xyz_file):

        with open(xyz_file, "w") as f:

            f.write()
        
            for at_id, at in self.atoms.items():

                self.a

        pass
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #























