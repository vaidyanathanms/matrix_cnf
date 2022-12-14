import numpy as np
import math

"""
A general definition of the single residue, common for all polysaccharides
    4
   /
  3 - 1
  \   /
    2  

    
"""



def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def pickAngle(index):
    if (index+1) % 4 == 1:
        angle = 0
    elif (index+1) % 4 == 2:
        angle = -1.5707963268
    elif (index+1) % 4 == 3:
        angle = -3.1415926536
    elif (index+1) % 4 == 0:
        angle = -4.7123889804
    return angle


class Polysaccharide:
    """
    The class that stores the parameters of a molecule composed of multiple residues.
    The class also has write_gro and write_top methods to create gro and top files
    """

    def __init__(self,name, unit_at):

        self.name      = name
        
        self.atoms     = []
        self.coordsX   = []
        self.coordsY   = []
        self.coordsZ   = []
        self.bonds     = []
        self.angles    = []
        self.dihedrals = []
        self.impropers = []
        
        self.qtot      = 0
        self.unit_at   = unit_at

    def gets_atoms(self, index, row):
        self.qtot += row[4]
        self.atoms.append([row[0]+index*self.unit_at, row[1], index+1, row[2], row[3], row[0]+index*self.unit_at, row[4], row[5], self.qtot])

    def gets_bonds(self, index, row):
        self.bonds.append([row[0]+index*self.unit_at,row[1]+index*self.unit_at,row[2],row[3]])

    def gets_angles(self, index, row):
        self.angles.append([row[0]+index*self.unit_at, row[1]+index*self.unit_at, row[2]+index*self.unit_at, row[3], row[4], row[5]])

    def gets_dihedrals(self, index, row):
        self.dihedrals.append([row[0]+index*self.unit_at, row[1]+index*self.unit_at, row[2]+index*self.unit_at, row[3]+index*self.unit_at, row[4], row[5], row[6], row[7]])
    
    def gets_impropers(self, index, row):
        self.impropers.append([row[0]+index*self.unit_at, row[1]+index*self.unit_at, row[2]+index*self.unit_at, row[3]+index*self.unit_at, row[4], row[5], row[6]])

    
    
    def writeGro(self):
        
        fname = self.name + ".gro"
        groFile = open(fname, "w")
        groFile.write("%s\n" % self.name)
        groFile.write(" %i\n" % len(self.atoms))
        count = 1

        for at, x, y, z  in zip(self.atoms,self.coordsX, self.coordsY, self.coordsZ):
            line = "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (at[2], at[3], at[4], count, x, y, z)
            count += 1
            if count == 100000:
                count = 0
            groFile.write(line)
        
        box = max([max(self.coordsX),max(self.coordsY),max(self.coordsZ)])
        boxX=boxY=boxZ=box+(box*0.2) 
        text = "%11.5f %11.5f %11.5f\n" % (boxX, boxY, boxZ )
        

        groFile.write(text)
    
        groFile.close()
    


    def writeTop(self):

        topText = []
        molName = self.name


        headTop = """
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
"""

        headMoleculetype = """
[ moleculetype ]
;name            nrexcl
%-16s 1
"""

        headAtomTypes = """
[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
"""

        headBonds = """
[ bonds ]
;   ai     aj funct   r             k
"""

        headAngles = """
[ angles ]
;   ai     aj     ak    funct   theta         cth
"""
        headProDih = """
[ dihedrals ] 
;    i      j      k      l   func    C0         C1         C2         C3         C4         C5
"""
        headImp = """
; impropers  
"""

        topText.append(headTop)
        topText.append(headMoleculetype % molName)
        topText.append(headAtomTypes)
        for row in self.atoms:
            line = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot %1.3f\n" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])
            topText.append(line)

        topText.append(headBonds)
        for row in self.bonds:
            line = "%6i %6i %3i %13.4e %13.4e\n" % (row[0], row[1], 1, row[2], row[3])
            topText.append(line)

        topText.append(headAngles)
        for row in self.angles:
            line = "%6i %6i %6i %6i %13.4e %13.4e\n" % (row[0], row[1], row[2], row[3], row[4], row[5])
            topText.append(line)

        topText.append(headProDih)
        for row in self.dihedrals:
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i\n" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])
            topText.append(line)

        topText.append(headImp)
        for row in self.impropers:
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f\n" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
            topText.append(line)

        footerTop = """
[ system ]
; Name
UNNAMED

[ molecules ]
; Compound        #mols
%-16s             1
"""

        topText.append(footerTop % molName)

        fileName = molName + ".top"
        topFile = open(fileName, "w")
        topFile.writelines(topText)
        topFile.close()

    
def GLCA14(units_num):

    name = "GLCA14"
    at = 4
    

    atom = []
    atom.append([1, "TP3", "AMY", "T1", 0, 36])
    atom.append([2, "TN4", "AMY", "T2", 0, 36])
    atom.append([3, "P3", "AMY", "R3", 0, 72])
    atom.append([4, "SN4", "AMY", "S4", 0, 54])
    
    coord = []

    coord.append([0.000, -0.2440, -0.0155])
    coord.append([0.000, -0.0040,  0.1045])
    coord.append([0.000,  0.2525,  0.0044])
    coord.append([0.000, -0.0040, -0.0935])

                
    bondIn = []
    bondIn.append([1, 4, 0.251, 18000])
    bondIn.append([2, 3, 0.280, 32000])
    bondIn.append([2, 4, 0.292, 60000])
    bondIn.append([3, 4, 0.279, 34000])

    bondConnect = []
    bondConnect.append([2, 8, 0.280, 32000])

    angleIn = []
    angleIn.append([1, 4, 2, 2,  86, 340])
    angleIn.append([1, 4, 3, 2, 142, 580])
    
    angleConnect = []
    angleConnect.append([3, 2, 8,  2, 103, 310]) 
    angleConnect.append([4, 2, 8, 10, 103, 180]) 
    angleConnect.append([1, 2, 8,  2, 105, 210]) 
    angleConnect.append([2, 8, 5,  2, 107,  30]) 
    angleConnect.append([2, 8, 6, 10, 148, 240]) 
    angleConnect.append([2, 8, 7,  2,  99, 120]) 
    
    dihIn = []
    
    
    dihConnect = []
    dihConnect.append([4, 2, 8, 6, 1, 161, -12.8, 1])
    dihConnect.append([3, 2, 8, 7, 1, 165, -2.2, 1])
    dihConnect.append([1, 3, 7, 5, 1, -72, -3.6, 5])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 9, 200])
    
    improperConnect = []
    improperConnect.append([8, 2, 7, 6, 2, 22.3, 170])
    improperConnect.append([2, 8, 4, 3, 2, -67, 300])



    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54

    axis = [1,0,0]
    for i in range(units_num):
        if (i+1) % 2 == 1:
            theta = 0
        elif (i+1) % 2 == 0:
            theta = 3.1415926536 # 180
        for row in coord:
            #row = rotation_matrix(axis, theta) @ row
            row = np.dot(rotation_matrix(axis, theta), row)
            mol.coordsX.append(row[0]+i*0.4)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
                mol.gets_impropers(0, row)

    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()
        

def GLCB12(units_num):

    name = "GLCB12"
    at = 4
    dist = 0.55

    atom = []
    atom.append([1, "TP3", "GLC2", "T1", 0, 36])
    atom.append([2, "TN4", "GLC2", "T2", 0, 36])
    atom.append([3, "P3",  "GLC2", "R3", 0, 72])
    atom.append([4, "SN4", "GLC2", "S4", 0, 54])

    coord = []

    coord.append([0.454, 0.200, 0.000])
    coord.append([0.275, 0.401, 0.000])
    coord.append([0.000, 0.400, 0.000])
    coord.append([0.202, 0.215, 0.000])

    bondIn = []
    bondIn.append([1, 4, 0.265, 10000])
    bondIn.append([2, 3, 0.263, 50000])
    bondIn.append([2, 4, 0.290, 50000])
    bondIn.append([3, 4, 0.272, 24000])

    bondConnect = []
    bondConnect.append([2, 7, 0.274, 6000])
    bondConnect.append([4, 7, 0.530, 600])


    angleIn = []
    angleIn.append([1, 4, 2, 2, 80, 200])
    angleIn.append([1, 4, 3, 2, 132, 400])
    

    angleConnect = []
    angleConnect.append([3, 2, 7, 2, 88, 30])
    angleConnect.append([1, 2, 7, 10, 165, 60])
    angleConnect.append([2, 7, 6, 2, 97, 50])
    angleConnect.append([2, 7, 8, 10, 148, 40]) 


    dihIn = []

    dihConnect = []
    dihConnect.append([3, 2, 7, 6, 1, -141, -22, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 11, 150])

    improperConnect = []
    improperConnect.append([7, 6, 8, 2, 2, 15, 300])
    improperConnect.append([2, 4, 3, 7, 2, 1, 200])

    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54
    
    for i in range(units_num):
        for row in coord:
            mol.coordsX.append(row[0]+i*dist)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()
        


def GLCA16(units_num):

    name = "GLCA16"
    at = 4 
    

    atom = []
    atom.append([1, "TP3", "GLC6", "T1", 0, 36])
    atom.append([2, "TN4", "GLC6", "T2", 0, 36])
    atom.append([3, "P3", "GLC6", "R3", 0, 72])
    atom.append([4, "SN4", "GLC6", "S4", 0, 54])

    coord = []

    coord.append([0.000, -0.2440, -0.0155])
    coord.append([0.000, -0.0040,  0.1045])
    coord.append([0.000,  0.2525,  0.0044])
    coord.append([0.000, -0.0040, -0.0935])


    bondIn = []
    bondIn.append([1, 4, 0.251, 21000])
    bondIn.append([2, 3, 0.280, 35000])
    bondIn.append([2, 4, 0.322, 60000])
    bondIn.append([3, 4, 0.289, 24000])

    bondConnect = []
    bondConnect.append([2, 5, 0.225, 13000])
    bondConnect.append([3, 7, 0.82, 2400])

    bondConnect3 = []
    bondConnect3.append([1, 9, 0.79, 180])
    
    angleIn = []
    angleIn.append([1, 4, 2, 2, 80, 600])
    angleIn.append([1, 4, 3, 2, 135, 560])

    angleConnect = []
    angleConnect.append([3, 2, 5, 2, 109, 180])
    angleConnect.append([4, 2, 5, 2, 105, 100])
    angleConnect.append([2, 5, 6, 10, 122, 28])

       

    dihIn = []
    
    dihConnect = []
    dihConnect.append([3, 2, 5, 6, 1, 155, -3, 2])
    dihConnect.append([2, 5, 6, 7, 1, 170, -3, 5])
    dihConnect.append([2, 5, 6, 7, 1, 60, 3, 1])
    dihConnect.append([2, 4, 8, 6, 1, -25, 5, 1])
    dihConnect.append([2, 4, 8, 6, 1, 165, -4, 3])

    improperIn = []
    improperIn.append([4, 2, 3, 1, 2, -7, 250])

    improperConnect = []
    improperConnect.append([2, 3, 5, 4, 2, -74, 110])


    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54

    axis = [1,0,0]
    for i in range(units_num):
        for row in coord:
            theta = pickAngle(i)
            row = np.dot(rotation_matrix(axis, theta), row)
            mol.coordsX.append(row[0]+i*0.4)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    if bondConnect3:
        for row in bondConnect3:
            for i in range(units_num - 2):
                mol.gets_bonds(i, row)


    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()



def GLCB13(units_num):

    name = "GLCB13"
    at = 4
    dist = 0.55

    atom = []
    atom.append([1, "TP3", "CUR", "T1", 0, 36])
    atom.append([2, "TN4", "CUR", "T2", 0, 36])
    atom.append([3, "P3", "CUR", "R3", 0, 72])
    atom.append([4, "SN4", "CUR", "S4", 0, 54])

    coord = []

    coord.append([0.454, 0.200, 0.000])
    coord.append([0.275, 0.401, 0.000])
    coord.append([0.000, 0.400, 0.000])
    coord.append([0.202, 0.215, 0.000])

    bondIn = []
    bondIn.append([1, 4, 0.268, 10000])
    bondIn.append([2, 3, 0.240, 38000])
    bondIn.append([2, 4, 0.287, 50000])
    bondIn.append([3, 4, 0.287, 24000])

    bondConnect = []
    bondConnect.append([2, 7, 0.250, 12000])


    angleIn = []
    angleIn.append([1, 4, 2, 2, 78, 220])
    angleIn.append([1, 4, 3, 2, 123, 400])

    angleConnect = []
    angleConnect.append([3, 2, 7, 2, 75, 60]) 
    angleConnect.append([1, 2, 7, 10, 148, 80]) 
    angleConnect.append([2, 7, 6, 2, 128, 110]) 
    angleConnect.append([2, 7, 8, 2, 65, 40]) 

    dihIn = []

    dihConnect = []
    dihConnect.append([3, 2, 7, 6, 1, -174, -20, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 15, 120])

    improperConnect = []
    improperConnect.append([7, 2, 8, 6, 2, 22, 200])
    improperConnect.append([2, 7, 3, 4, 2, 3, 200])

    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54

    
    for i in range(units_num):
        for row in coord:
            mol.coordsX.append(row[0]+i*dist)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 
    

    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
            for i in range(units_num):
                mol.gets_impropers(i, row)
    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()


def GLCB14(units_num):

    name = "GLCB14"
    at = 4
    

    atom = []
    atom.append([1, "TP3", "CELL", "T1", 0, 36])
    atom.append([2, "TN4", "CELL", "T2", 0, 36])
    atom.append([3, "P3", "CELL", "R3", 0, 72])
    atom.append([4, "SN4", "CELL", "S4", 0, 54])
    
    coord = []

    coord.append([0.000, -0.2440, -0.0155])
    coord.append([0.000, -0.0040,  0.1045])
    coord.append([0.000,  0.2525,  0.0044])
    coord.append([0.000, -0.0040, -0.0935])

    

    bondIn = []
    bondIn.append([1, 4, 0.250, 14100])
    bondIn.append([2, 3, 0.2680, 37500])
    bondIn.append([2, 4, 0.2570, 53200])
    bondIn.append([3, 4, 0.2730, 27000])

    bondConnect = []
    bondConnect.append([2, 8, 0.267, 7500])
    bondConnect.append([2, 6, 0.520, 16300])
    bondConnect.append([4, 8, 0.542, 3770])

    angleIn = []
    angleIn.append([1, 4, 2, 2, 91, 220])
    angleIn.append([1, 4, 3, 10, 143, 159])
    
    angleConnect = []
    angleConnect.append([3, 2, 8, 2, 115, 245]) 
    angleConnect.append([1, 2, 8, 2, 127, 350]) 
    angleConnect.append([2, 8, 5, 2, 123, 16]) 
    angleConnect.append([2, 8, 7, 2, 93, 52]) 
    
    dihIn = []
    

    dihConnect = []
    dihConnect.append([3, 2, 8, 7, 1, -135 ,-35, 1])

    improperIn = []
    improperIn.append([4, 3, 2, 1, 2, 9, 200])

    improperConnect = []
    improperConnect.append([2, 3, 8, 4, 2, 9, 229])
    improperConnect.append([8, 2, 7, 6, 2, 11, 212])   


    mol = Polysaccharide(name, at)

    for i in range(units_num):
            for row in atom:
                mol.gets_atoms(i, row)

    mol.atoms[-3][1] = "SN4"
    mol.atoms[-3][4] = "S2"
    mol.atoms[-3][7] = 54
        

    axis = [1,0,0]
    for i in range(units_num):
        if (i+1) % 2 == 1:
            theta = 0
        elif (i+1) % 2 == 0:
            theta = 3.1415926536 # 180
        for row in coord:
            row = np.dot(rotation_matrix(axis, theta), row)
            mol.coordsX.append(row[0]+i*0.4)
            mol.coordsY.append(row[1]) 
            mol.coordsZ.append(row[2]) 


    if bondIn:
        for row in bondIn:
            for i in range(units_num):
                mol.gets_bonds(i, row)

    if bondConnect:
        for row in bondConnect:
            for i in range(units_num - 1):
                mol.gets_bonds(i, row)

    
    if angleIn:
        for row in angleIn:
            for i in range(units_num):
                mol.gets_angles(i, row)
    if angleConnect:
        for row in angleConnect:
            for i in range(units_num-1):
                mol.gets_angles(i, row)

    if dihIn:
        for row in dihIn:
            for i in range(units_num):
                mol.gets_dihedrals(i, row)
    if dihConnect:
        for row in dihConnect:
            for i in range(units_num-1):
                mol.gets_dihedrals(i, row)

        
    if improperIn:
        for row in improperIn:
                mol.gets_impropers(0, row)

    if improperConnect:
        for row in improperConnect:
            for i in range(units_num-1):
                mol.gets_impropers(i, row)

    
    mol.writeGro()
    mol.writeTop()


def main():

    units_name = ["alpha(1-4)", "beta(1-2)", "alpha(1-6)", "beta(1-4)", "beta(1-3)"]
    print("Select linkage type:")
    for i , name in enumerate(units_name, 1):
        print("%d %s" % (i, name))
    sugar = int(input())
    print("You chose %s" % (units_name[sugar - 1]))
    while True:
        print("Insert number of residue in a homopolysaccharide:")
        units_nr = int(input())
        if(units_nr==1):
            print("The number of residues must be larger than 1")
            continue
        else:
            break
    if(sugar == 1):
        GLCA14(units_nr)
    if(sugar == 2):
        GLCB12(units_nr)
    if(sugar == 3):
        GLCA16(units_nr)
    if(sugar == 4):
        GLCB14(units_nr)
    if(sugar == 5):
        GLCB13(units_nr)
    

if __name__ == "__main__":
    main()
