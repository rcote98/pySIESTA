
import numpy as np

from pySIESTA.constants import ang2bohr, e2C, bohr2m
from pySIESTA.structures import SolidGeometry
from pySIESTA.scheduler import FDFSettings

class STOGeometry(SolidGeometry):

    def __init__(self, supercell, lattice=3.874,  reference="cubic"):

        species = {
            0: [1, 38, "Sr"], 
            1: [2, 22, "Ti"], 
            2: [3, 8,  "Ox"],
            3: [3, 8,  "Oy"],
            4: [3, 8,  "Oz"],  
        }

        super().__init__(supercell, species)
        self.loaded_struct = True

        for i in range(3):
            self.lat_vecs[i,i] = lattice*self.supercell[i]
            self.ref_lat_vecs[i,i] = lattice*self.supercell[i]

        self.born_charges = {
            0: np.array([ 2.5256,  2.5256,  2.5256]), 
            1: np.array([ 7.5526,  7.5526,  7.5526]), 
            2: np.array([-5.9494, -2.0644, -2.0644]),
            3: np.array([-2.0644, -5.9494, -2.0644]),
            4: np.array([-2.0644, -2.0644, -5.9494]) 
        }

        sc = self.supercell

        uc_pos = {
            0: np.array([0.00, 0.00, 0.00]),
            1: np.array([0.50, 0.50, 0.50]),
            2: np.array([0.00, 0.50, 0.50]),
            3: np.array([0.50, 0.00, 0.50]),
            4: np.array([0.50, 0.50, 0.00])
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
        Ox, Oy, Oz = 2, 3, 4

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

        A, B, Ox, Oy, Oz = 0, 1, 2, 3, 4

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class STO_FDFSettings(FDFSettings):

    
    def __init__(self):

        # same as the normal fdf class but with preloaded stuff

        self.settings = {
            'systemname': ['Strontium titanate (SrTiO3)'], 
            'systemlabel': ['srtio3-str'], 
            'numberofspecies': [3], 
            'numberofatoms': [40], 
            'chemicalspecieslabel': [
                    [1, 38, 'Sr'], 
                    [2, 22, 'Ti'], 
                    [3, 8, 'O']
            ], 
            'ps.lmax': [
                    ['Sr', 3], 
                    ['Ti', 3], 
                    ['O',  3]], 
            'pao.basis': [
                    ['sr', '5', '1.64000'], 
                    ['n=4', '0', '1', 'e', '155.00000', '6.00000'], 
                    ['6.49993'], 
                    ['1.00000'], 
                    ['n=5', '0', '2', 'e', '149.48000', '6.50000'], 
                    ['6.99958', '5.49957'], 
                    ['1.00000', '1.00000'], 
                    ['n=4', '1', '1', 'e', '148.98000', '5.61000'], 
                    ['6.74964'], 
                    ['1.00000'], 
                    ['n=5', '1', '1', 'e', '4.57000', '1.20000'], 
                    ['4.00000'], 
                    ['1.00000'], 
                    ['n=4', '2', '1', 'e', '146.26000', '6.09000'], 
                    ['6.63062'], 
                    ['1.00000'], 
                    ['ti', '5', '1.91'], 
                    ['n=3', '0', '1', 'e', '93.95', '5.20'], 
                    ['5.69946662616249'], 
                    ['1.00000000000000'], 
                    ['n=3', '1', '1', 'e', '95.47', '5.20'], 
                    ['5.69941339465994'], 
                    ['1.00000000000000'], 
                    ['n=4', '0', '2', 'e', '96.47', '5.60'], 
                    ['6.09996398975307', '5.09944363262274'], 
                    ['1.00000000000000', '1.00000000000000'], 
                    ['n=3', '2', '2', 'e', '46.05', '4.95'], 
                    ['5.94327035784617', '4.70009988294302'], 
                    ['1.00000000000000', '1.00000000000000'], 
                    ['n=4', '1', '1', 'e', '0.50', '1.77'], 
                    ['3.05365979938936'], 
                    ['1.00000000000000'], 
                    ['o', '3', '-0.28'], 
                    ['n=2', '0', '2', 'e', '40.58', '3.95'], 
                    ['4.95272270428712', '3.60331408800389'], 
                    ['1.00000000000000', '1.00000000000000'], 
                    ['n=2', '1', '2', 'e', '36.78', '4.35'], 
                    ['4.99990228025066', '3.89745395068600'], 
                    ['1.00000000000000', '1.00000000000000'], 
                    ['n=3', '2', '1', 'e', '21.69', '0.93'], 
                    ['2.73276990670788'], 
                    ['1.00000000000000']
                    ], 
            'writecoorstep': ['.true.'], 
            'AtomicCoordinatesFormat': ['fractional'],
            'AtomCoorFormatOut': ['fractional'],
            'kgrid_monkhorst_pack': [
                    ['3', '0', '0', '0.5'], 
                    ['0', '3', '0', '0.5'], 
                    ['0', '0', '3', '0.5']
                    ], 
            'xc.functional': ['lda'], 
            'xc.authors': ['ca'], 
            'meshcutoff': [600, 'ry'], 
            'dm.numberpulay': [3], 
            'dm.usesavedm': ['.true.'], 
            'dm.tolerance': ['1.d-4'], 
            'maxscfiterations': ['100'], 
            'electronictemperature': ['0.025', 'ev'], 
            'scf.mixafterconvergence': ['.false.'], 
            'md.typeofrun': ['cg'], 
            'md.variablecell': ['.true.'], 
            'md.usestructfile': ['.true.'], 
            'md.usesavexv': ['.false.'], 
            'md.usesavecg': ['.false.'], 
            'md.numcgsteps': [40], 
            'md.maxcgdispl': [0.3, 'bohr'], 
            'md.maxforcetol': [0.01, 'ev/ang'], 
            'md.maxstresstol': [0.0001, 'ev/ang**3'], 
            'geometryconstraints': [
                ['routine', 'constr']
                ], 
            'projecteddensityofstates': [
                ['-70.00', '5.00', '0.050', '3000', 'ev']
                ], 
            'pdos.kgrid_monkhorst_pack': [
                    ['30', '0', '0', '0.5'], 
                    ['0', '30', '0', '0.5'], 
                    ['0', '0', '30', '0.5']
                    ]
        }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def STO_AFD(supercell, angle, mode="a", axis="z", clockwise=False):

    """

    """

    _, _, Ox, Oy, Oz = 0, 1, 2, 3, 4
    sc = np.array(supercell)
    geom = STOGeometry(sc)
    disp = 0.5*np.arctan(angle/180.*np.pi)/sc

    if clockwise:
        cw = 1
    else:
        cw = 0

    for x in range(supercell[0]):
        for y in range(supercell[1]):
            for z in range(supercell[2]):

                if mode == "a":

                    factor = (-1)**x * (-1)**y * (-1)**z * (-1)**cw

                    if axis == "x":
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                    elif axis == "y":
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                    elif axis == "z":
                        geom.positions[x,y,z,Ox,1] -= factor*disp[1]
                        geom.positions[x,y,z,Oy,0] += factor*disp[0]
                    elif axis == "xy" or axis == "yx":
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                    else:
                        raise NotImplementedError()

                elif mode == "i":

                    if axis == "x":
                        factor = (-1)**y * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                    elif axis == "y":
                        factor = (-1)**x * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                    elif axis == "z":
                        factor = (-1)**x * (-1)**y * (-1)**cw
                        geom.positions[x,y,z,Ox,1] -= factor*disp[1]
                        geom.positions[x,y,z,Oy,0] += factor*disp[0]
                    elif axis == "xy" or axis == "yx":
                        factor = (-1)**y * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                        factor = (-1)**x * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                    else:
                        raise NotImplementedError()
    
                else:
                    raise NotImplementedError()

    return geom


def STO_FE(supercell, disp, axis="z"):

    """

    """

    _, B, _, _, _ = 0, 1, 2, 3, 4
    sc = np.array(supercell)
    geom = STOGeometry(sc)


    for x in range(supercell[0]):
        for y in range(supercell[1]):
            for z in range(supercell[2]):

                if axis == "x":
                    geom.positions[x,y,z,B,0] += disp
                elif axis == "y":
                    geom.positions[x,y,z,B,1] += disp
                elif axis == "z":
                    geom.positions[x,y,z,B,2] += disp
                elif axis == "xy" or axis == "yx":
                    geom.positions[x,y,z,B,0] += disp
                    geom.positions[x,y,z,B,1] += disp
                else:
                    raise NotImplementedError()

    return geom


def STO_AFD_FE(supercell, angle, ti_disp, mode="a", axis="z", clockwise=False):

    """

    """

    _, B, Ox, Oy, Oz = 0, 1, 2, 3, 4
    sc = np.array(supercell)
    geom = STOGeometry(sc)
    disp = 0.5*np.arctan(angle/180.*np.pi)/sc

    if clockwise:
        cw = 1
    else:
        cw = 0

    for x in range(supercell[0]):
        for y in range(supercell[1]):
            for z in range(supercell[2]):

                if mode == "a":

                    factor = (-1)**x * (-1)**y * (-1)**z * (-1)**cw

                    if axis == "x":
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                        geom.positions[x,y,z,B,0]  += ti_disp
                    elif axis == "y":
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                        geom.positions[x,y,z,B,1]  += ti_disp
                    elif axis == "z":
                        geom.positions[x,y,z,Ox,1] -= factor*disp[1]
                        geom.positions[x,y,z,Oy,0] += factor*disp[0]
                        geom.positions[x,y,z,B,2]  += ti_disp
                    elif axis == "xy" or axis == "yx":
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                        geom.positions[x,y,z,B,0]  += ti_disp
                        geom.positions[x,y,z,B,1]  += ti_disp
                    else:
                        raise NotImplementedError()

                elif mode == "i":

                    if axis == "x":
                        factor = (-1)**y * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                        geom.positions[x,y,z,B,0]  += ti_disp
                    elif axis == "y":
                        factor = (-1)**x * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                        geom.positions[x,y,z,B,1]  += ti_disp
                    elif axis == "z":
                        factor = (-1)**x * (-1)**y * (-1)**cw
                        geom.positions[x,y,z,Ox,1] -= factor*disp[1]
                        geom.positions[x,y,z,Oy,0] += factor*disp[0]
                        geom.positions[x,y,z,B,2]  += ti_disp
                    elif axis == "xy" or axis == "yx":
                        factor = (-1)**y * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oy,2] -= factor*disp[2]
                        geom.positions[x,y,z,Oz,1] += factor*disp[1]
                        factor = (-1)**x * (-1)**z * (-1)**cw
                        geom.positions[x,y,z,Oz,0] -= factor*disp[0]
                        geom.positions[x,y,z,Ox,2] += factor*disp[2]
                        geom.positions[x,y,z,B,0]  += ti_disp
                        geom.positions[x,y,z,B,1]  += ti_disp
                    else:
                        raise NotImplementedError()
    
                else:
                    raise NotImplementedError()

    return geom