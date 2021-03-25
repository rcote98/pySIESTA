
import numpy as np

from pySIESTA.structures import SolidGeometry
from pySIESTA.constants import ang2bohr, e2C, bohr2m

class STOGeometry(SolidGeometry):

    def __init__(self, supercell, lattice=3.874,  reference="cubic"):

        species = {
            0: [1, 38, "Sr"], 
            1: [2, 22, "Ti"], 
            4: [3, 8,  "Ox"],
            3: [3, 8,  "Oy"],
            2: [3, 8,  "Oz"],  
        }

        super().__init__(supercell, species)
        self.loaded_struct = True

        for i in range(3):
            self.lat_vecs[i,i] = lattice*self.supercell[i]
            self.ref_lat_vecs[i,i] = lattice*self.supercell[i]

        self.born_charges = {
            0: np.array([ 2.5256,  2.5256,  2.5256]), 
            1: np.array([ 7.5526,  7.5526,  7.5526]), 
            4: np.array([-5.9494, -2.0644, -2.0644]),
            3: np.array([-2.0644, -5.9494, -2.0644]),
            2: np.array([-2.0644, -2.0644, -5.9494]) 
        }

        sc = self.supercell

        uc_pos = {
            0: np.array([0.00, 0.00, 0.00]),
            1: np.array([0.50, 0.50, 0.50]),
            4: np.array([0.00, 0.50, 0.50]),
            3: np.array([0.50, 0.00, 0.50]),
            2: np.array([0.50, 0.50, 0.00])
        }

        for x in range(self.supercell[0]):
            for y in range(self.supercell[1]):
                for z in range(self.supercell[2]):
                    for at in range(self.nats):
                        self.positions[x,y,z,at,0] = uc_pos[at][0]/sc[0] + x/sc[0]
                        self.positions[x,y,z,at,1] = uc_pos[at][1]/sc[1] + y/sc[1]
                        self.positions[x,y,z,at,2] = uc_pos[at][2]/sc[2] + z/sc[2]

        self.reference = np.copy(self.positions)
        self.loaded_ref = True

    def rotations(self, angles = True):

        sc = self.supercell
        Ox, Oy, Oz = 4, 3, 2

        ROT_X=[
            # atom, hopping, weight, target vector
            [Oz, [ 0, 0, 0],  1/2., [ 0.0, 0.7071067812, 0.0]],
            [Oy, [ 0, 0, 0], -1/2., [ 0.0, 0.0, 0.7071067812]],
            [Oz, [ 0, 0, 1], -1/2., [ 0.0, 0.7071067812, 0.0]],
            [Oy, [ 0, 1, 0],  1/2., [ 0.0, 0.0, 0.7071067812]],
        ]

        ROT_Y=[
            # atom, hopping, weight, target vector
            [Ox, [0, 0, 0], 1/2.,[ 0.0, 0.0, 0.7071067812]],
            [Oz, [0, 0, 0],-1/2.,[ 0.7071067812, 0.0, 0.0]],
            [Ox, [1, 0, 0],-1/2.,[ 0.0, 0.0, 0.7071067812]],
            [Oz, [0, 0, 1], 1/2.,[ 0.7071067812, 0.0, 0.0]],
        ]

        ROT_Z=[
            # atom, hopping, weight, target vector
            [Ox, [ 0, 0, 0], -1/2., [ 0.0, 0.7071067812, 0.0]],
            [Oy, [ 0, 0, 0],  1/2., [ 0.7071067812, 0.0, 0.0]],
            [Ox, [ 1, 0, 0],  1/2., [ 0.0, 0.7071067812, 0.0]],
            [Oy, [ 0, 1, 0], -1/2., [ 0.7071067812, 0.0, 0.0]],
        ]
        
        rots = np.zeros((sc[0], sc[1], sc[2], 3))
        disps = self.get_displacements(strain = False)

        for x in range(sc[0]):
            for y in range(sc[1]):
                for z in range(sc[2]):

                    cell = np.array([x,y,z])
                
                    for atom in ROT_X:
                        atom_cell = np.mod(cell + atom[1], sc)
                        nx, ny, nz = atom_cell
                        rots[x,y,z,0] += atom[2]*np.dot(atom[3], disps[nx,ny,nz,atom[0],:])

                    for atom in ROT_Y:
                        atom_cell = np.mod(cell + atom[1], sc)
                        nx, ny, nz = atom_cell
                        rots[x,y,z,1] += atom[2]*np.dot(atom[3], disps[nx,ny,nz,atom[0],:])

                    for atom in ROT_Z:
                        atom_cell = np.mod(cell + atom[1], sc)
                        nx, ny, nz = atom_cell
                        rots[x,y,z,2] += atom[2]*np.dot(atom[3], disps[nx,ny,nz,atom[0],:])

        if angles:
            BOdist = self.ref_lat_vecs[0,0]/sc[0]/2
            return np.arctan(rots/(np.sqrt(2)*BOdist))*180/np.pi
        else:
            return rots


    def polarization(self):

        A, B, Ox, Oy, Oz = 0, 1, 4, 3, 2

        FE_mode=[  # atom, hopping, weight
            # "frame"
            [A, [0, 0, 0], 1./8.],
            [A, [1, 0, 0], 1./8.],
            [A, [1, 1, 0], 1./8.],
            [A, [0, 1, 0], 1./8.],
            [A, [0, 0, 1], 1./8.],
            [A, [1, 0, 1], 1./8.],
            [A, [1, 1, 1], 1./8.],
            [A, [0, 1, 1], 1./8.],
            # "octahedra"
            [B, [0, 0, 0], 1.], # b site
            [Ox, [0, 0, 0], 1./2.], 
            [Ox, [1, 0, 0], 1./2.],
            [Oy, [0, 0, 0], 1./2.],
            [Oy, [0, 1, 0], 1./2.],
            [Oz, [0, 0, 0], 1./2.],
            [Oz, [0, 0, 1], 1./2.]
        ]

        sc = self.supercell
        uc_vol = np.prod(np.linalg.norm(self.lat_vecs*ang2bohr, axis=0))/np.prod(self.supercell)

        pols = np.zeros((sc[0], sc[1], sc[2], 3))
        disps = self.get_displacements(unit="bohr", strain = False)

        for x in range(sc[0]):
            for y in range(sc[1]):
                for z in range(sc[2]):

                    cell = np.array([x,y,z])

                    for atom in FE_mode:

                        atom_cell = np.mod(cell + atom[1], sc)
                        nx, ny, nz = atom_cell

                        charge = np.array(self.born_charges[atom[0]])

                        pols[x,y,z,0] += atom[2]*charge[0]*disps[nx,ny,nz,atom[0],0]
                        pols[x,y,z,1] += atom[2]*charge[1]*disps[nx,ny,nz,atom[0],1]
                        pols[x,y,z,2] += atom[2]*charge[2]*disps[nx,ny,nz,atom[0],2]

        pols = pols/uc_vol          # in e/bohr2
        pols = pols*e2C/bohr2m**2   # in C/m2

        return pols
