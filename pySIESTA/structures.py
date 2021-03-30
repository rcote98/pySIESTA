"""
    Provides data structures to handle SIESTA files.
"""


import csv

# third party imports
import numpy as np

from pySIESTA.constants import ang2bohr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODULE STRUCTURE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
# + class SolidGeometry()
#   - reset()
#   - read_STRUCT(struct_file)
#   - write_XYZ(xyz_file)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class SolidGeometry():
   
    """
    SIESTA Geometry Container

    Attributes:
    ----------

     - atoms (dict): dictionary with positions of the atoms (fractional)
     - lat_vecs (array): supercell shape 

    """

    def __init__(self, supercell, species):

        self.species    = species
        self.nats       = len(species)
        self.supercell  = np.array(supercell)
        self.tot_ats    = np.prod(self.supercell)*self.nats 
        self.reset()

    def reset(self):

        sc = self.supercell
        self.loaded_struct  = False
        self.loaded_ref     = False
        self.strain         = np.identity(3) 
        self.lat_vecs       = np.zeros((3,3))
        self.ref_lat_vecs   = np.zeros((3,3))
        self.positions      = np.zeros((sc[0], sc[1], sc[2], self.nats, 3))
        self.reference      = np.zeros((sc[0], sc[1], sc[2], self.nats, 3))
        self.born_charges   = {}
        for at in self.species.keys():
            self.born_charges[at] = np.array([0.0, 0.0, 0.0])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def add_reference(self, struct_file):

        with open(struct_file, "r") as f:
            
            for i in range(3):
                line = f.readline().split()
                self.ref_lat_vecs[i,:] = [float(x) for x in line]

            self.strain = np.dot(self.lat_vecs, np.linalg.inv(self.ref_lat_vecs))

            if self.tot_ats != int(f.readline().split()[0]):
                print("ERROR: Different number of atoms in reference file.")
                exit

            for x in range(self.supercell[0]):
                for y in range(self.supercell[1]):
                    for z in range(self.supercell[2]):
                        for at in range(self.nats):

                            line = f.readline().split()
                            self.reference[x,y,z,at,:] = np.array([float(x) for x in line[2:]])

        self.loaded_ref = True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def set_strain(self, strain_mat):

        self.strain = strain_mat
        self.lat_vecs = np.dot(self.strain, self.ref_lat_vecs)
        
        pass

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def get_positions(self, unit = "angstrom"):

        sc = self.supercell
        tpos = np.zeros((sc[0], sc[1], sc[2], self.nats, 3))

        for x in range(self.supercell[0]):
            for y in range(self.supercell[1]):
                for z in range(self.supercell[2]):
                    for at in range(self.nats):
                        for d in range(3):
                            tpos[x,y,z,at,:] += self.positions[x,y,z,at,d]*self.lat_vecs[d,:]

        if unit == "angstrom":
            return tpos
        elif unit == "bohr":
            return tpos*ang2bohr
        else:
            print("ERROR: Unrecognized unit.")
            exit

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def get_displacements(self, strain = False, unit = "angstrom"):

        if not self.loaded_ref:
            print("WARNING: No reference structure loaded.")

        if strain: 

            sc = self.supercell
            tpos = np.zeros((sc[0], sc[1], sc[2], self.nats, 3))
            tref = np.zeros((sc[0], sc[1], sc[2], self.nats, 3))

            for x in range(self.supercell[0]):
                for y in range(self.supercell[1]):
                    for z in range(self.supercell[2]):
                        for at in range(self.nats):
                            for d in range(3):
                                tpos[x,y,z,at,:] += self.positions[x,y,z,at,d]*self.lat_vecs[d,:]
                                tref[x,y,z,at,:] += self.reference[x,y,z,at,d]*self.lat_vecs[d,:]

            tdisp = tpos - tref

        else:

            sc = self.supercell
            tdisp = np.zeros((sc[0], sc[1], sc[2], self.nats, 3))
            disps = self.positions - self.reference

            for x in range(self.supercell[0]):
                for y in range(self.supercell[1]):
                    for z in range(self.supercell[2]):
                        for at in range(self.nats):
                            for d in range(3):
                                tdisp[x,y,z,at,:] += disps[x,y,z,at,d]*self.ref_lat_vecs[d,:]

        if unit == "angstrom":
            return tdisp
        elif unit == "bohr":
            return tdisp*ang2bohr
        else:
            print("ERROR: Unrecognized unit.")
            exit

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def polarization(self):

        if not self.loaded_ref:
            print("WARNING: No reference structure loaded.")

        sc = self.supercell
        pols = np.zeros((sc[0], sc[1], sc[2], 3))

        return pols

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def read_STRUCT(self, struct_file):

        with open(struct_file, "r") as f:
            
            if not self.loaded_ref:

                for i in range(3):
                    line = f.readline().split()
                    self.lat_vecs[i,:] = [float(x) for x in line]
                    self.ref_lat_vecs[i,:] = [float(x) for x in line]

                if self.tot_ats != int(f.readline().split()[0]):
                    print("ERROR: Different number of atoms in structure file.")
                    exit

                for x in range(self.supercell[0]):
                    for y in range(self.supercell[1]):
                        for z in range(self.supercell[2]):
                            for at in range(self.nats):
                                line = f.readline().split()
                                self.positions[x,y,z,at,:] = np.array([float(x) for x in line[2:]])

            else:

                for i in range(3):
                    line = f.readline().split()
                    self.lat_vecs[i,:] = [float(x) for x in line]

                self.strain = np.dot(self.lat_vecs, np.linalg.inv(self.ref_lat_vecs))

                if self.tot_ats != int(f.readline().split()[0]):
                    print("ERROR: Different number of atoms in structure file.")
                    exit

                for x in range(self.supercell[0]):
                    for y in range(self.supercell[1]):
                        for z in range(self.supercell[2]):
                            for at in range(self.nats):
                                line = f.readline().split()
                                self.positions[x,y,z,at,:] = np.array([float(x) for x in line[2:]])

                pass

        self.loaded_struct = True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def write_STRUCT(self, struct_file):

        with open(struct_file, "w") as f:

            tsv = csv.writer(f, delimiter="\t")

            for i in range(3):
                line = ["{:10.8F}".format(d) for d in self.lat_vecs[i,:]] + ["{:10.8F}".format(0)]*3
                tsv.writerow(line)

            tsv.writerow([self.tot_ats])
            
            for x in range(self.supercell[0]):
                for y in range(self.supercell[1]):
                    for z in range(self.supercell[2]):
                        for at in range(self.nats):

                            line = []
                            line.append(self.species[at][0])
                            line.append(self.species[at][1])
                            line += ["{:10.8F}".format(d) for d in self.positions[x,y,z,at,:]]
                            tsv.writerow(line)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #






















