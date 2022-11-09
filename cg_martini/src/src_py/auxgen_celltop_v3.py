# To generate itp files for certain polymers and acetylated cellulose
# Using Martini V3.0 for cellulose structure
# Main file: gen_top.py
# Reference_1: Lutsyk et al., JCTC, 18, 5089-5017, 2022
# Table 4 of Reference 1 has the interaction parameters for beta(1-4)
# https://github.com/ricalessandri/Martini3-small-molecules/blob/main/models/martini_v3.0.0_small_molecules_v1.itp
# Reference_2: https://pubs.acs.org/doi/10.1021/acs.jctc.2c00553.
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
#----------------------------------------------------------------------------------
# Copied from carbo2martini.py (Uses a class to generate systems)
# gro file is generated from atomistic coordinates. So no need to worry about that
"""
A general definition of the single residue, common for all polysaccharides

    1
   /
  4 - 2
  \   /
    3  

    
"""
class Polysaccharide:
    """
    The class that stores the parameters of a molecule composed of multiple residues.
    The class also has write_gro and write_top methods to create gro and top files
    """

    def __init__(self,name, unit_at):

        self.name      = name
        self.atoms     = []
        self.bonds     = []
        self.angles    = []
        self.dihedrals = []
        self.impropers = []
        self.qtot      = 0
        self.unit_at   = unit_at

    def gets_atoms(self, index, row, prev_endID,flag1):
        self.qtot += row[4]
        self.atoms.append([row[0]+index*self.unit_at+prev_endID+flag1,\
                           row[1], index+1, row[2], row[3], \
                           row[0]+index*self.unit_at+prev_endID+flag1,\
                           row[4], row[5], self.qtot])

    def gets_bonds(self, index, row, prev_endID, flag1, flag2):
        self.bonds.append([row[0]+index*self.unit_at+prev_endID+flag1,\
                           row[1]+index*self.unit_at+prev_endID+flag2,\
                           row[2],row[3]])

    def gets_angles(self, index, row, prev_endID, flag1, flag2, flag3):
        self.angles.append([row[0]+index*self.unit_at+prev_endID+flag1,\
                            row[1]+index*self.unit_at+prev_endID+flag2,\
                            row[2]+index*self.unit_at+prev_endID+flag3,\
                            row[3], row[4], row[5]])

    def gets_dihedrals(self, index, row, prev_endID, flag1, flag2, flag3, flag4):
        self.dihedrals.append([row[0]+index*self.unit_at+prev_endID+flag1,\
                               row[1]+index*self.unit_at+prev_endID+flag2,\
                               row[2]+index*self.unit_at+prev_endID+flag3,\
                               row[3]+index*self.unit_at+prev_endID+flag4,\
                               row[4], row[5], row[6], row[7]])
    
    def gets_impropers(self, index, row, prev_endID, flag1, flag2, flag3, flag4):
        self.impropers.append([row[0]+index*self.unit_at+prev_endID+flag1,\
                               row[1]+index*self.unit_at+prev_endID+flag2,\
                               row[2]+index*self.unit_at+prev_endID+flag3,\
                               row[3]+index*self.unit_at+prev_endID+flag4,\
                               row[4], row[5], row[6]])

    def writeTop(self,fname,namemol): # define topology files

        topText = []
        molName = namemol #self.name

        headTop = """
; Generated using carbo2martini.py and gen_top.py
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
            line = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot %1.3f\n"\
                   % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])
            topText.append(line)

        topText.append(headBonds)
        for row in self.bonds:
            line = "%6i %6i %3i %13.4e %13.4e\n" \
                   % (row[0], row[1], 1, row[2], row[3])
            topText.append(line)

        topText.append(headAngles)
        for row in self.angles:
            line = "%6i %6i %6i %6i %13.4e %13.4e\n" \
                   % (row[0], row[1], row[2], row[3], row[4], row[5])
            topText.append(line)

        topText.append(headProDih)
        for row in self.dihedrals:
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i\n" \
                   % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])
            topText.append(line)

        topText.append(headImp)
        for row in self.impropers:
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f\n"\
                   % (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
            topText.append(line)

        footerTop = """
[ system ]
; Name
Combined
        
[ molecules ]
; Compound        #mols
%-16s             1
"""

        topText.append(footerTop % molName)

        fileName = fname + ".top"
        topFile = open(fileName, "w")
        topFile.writelines(topText)
        topFile.close()
#----------------------------------------------------------------------------------
# Define cellulose beta 1-4 version. Taken from carbo2martini3.py
def GLCB14(units_num,ncnf_bundles,ch_per_cnf,glycan_list,outfname,molname):

    name = "GLCB14"
    at = 4 # this defines the number of beads per monomer
    

    atom = []
    atom.append([1, "TP3", "CELL", "T1", 0, 36])
    atom.append([2, "TN4", "CELL", "T2", 0, 36])
    atom.append([3, "P3", "CELL", "R3", 0, 72])
    atom.append([4, "SN4", "CELL", "S4", 0, 54])
    atom.append([5, "SN4a", "CELL", "A1", 0, 54]) # Acetylation - Check params

    bondIn = []
    bondIn.append([1, 4, 0.250, 14100])
    bondIn.append([2, 3, 0.2680, 37500])
    bondIn.append([2, 4, 0.2570, 53200])
    bondIn.append([3, 4, 0.2730, 27000])
    bondIn.append([1, 5, 0.2730, 27000]) #check parameters

    bondConnect = []
    bondConnect.append([2, 8, 0.267, 7500])
    bondConnect.append([2, 6, 0.520, 16300])
    bondConnect.append([4, 8, 0.542, 3770])

    angleIn = []
    angleIn.append([1, 4, 2, 2, 91, 220])
    angleIn.append([1, 4, 3, 10, 143, 159])
    angleIn.append([5, 1, 4, 10, 143, 159]) #check parameters

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

    acet_flag   = []; acet_cnt = 0
    prev_chend_id = 0 #keeps track of number of beads in a given chain 
    for j in range(ncnf_bundles*ch_per_cnf):
        bead_cnt = 0
        for i in range(units_num): #units_num is same as nresidues or nmons
            glycan_index = i + (j-1)*units_num
            bead_cnt += int(0.5*len(glycan_list[glycan_index])) #0.5 for 2 ids per bead
            if i == 0 and j == 0: #flag for first monomer is zero
                if len(glycan_list[glycan_index]) == 8:
                    acet_flag.append(acet_cnt)
                    for k in range(0,4):
                        mol.gets_atoms(i,atom[k],prev_chend_id,0)  
                elif len(glycan_list[glycan_index]) == 10:
                    acet_cnt += 1
                    acet_flag.append(acet_cnt)
                    for k in range(0,5):
                        mol.gets_atoms(i,atom[k],prev_chend_id,0)
                elif len(glycan_list[glycan_index]) == 12:
                    acet_cnt += 2
                    acet_flag.append(acet_cnt)
                    for k in range(0,6):
                        mol.gets_atoms(i,atom[k],prev_chend_id,0)
                else:
                    print(glycan_list[glycan_index])
                    raise RuntimeError("More than 5 atoms in glycan for "\
                                       + str(glycan_index))
            else: # flag for other monomers will depend on acet_cnt in previous monomer
                # 4 beads per residue or unit or mon and 2 ids per bead
                if len(glycan_list[glycan_index]) == 8:
                    acet_flag.append(acet_cnt)
                    for k in range(0,4):
                        mol.gets_atoms(i,atom[k],prev_chend_id,acet_flag[i-1])
                    if i == units_num-1: #Changes to end monomers 
                        mol.atoms[-3][1] = "SN4"
                        mol.atoms[-3][4] = "S2"
                        mol.atoms[-3][7] = 54
                # 5 beads per residue or unit or mon and 2 ids per bead
                elif len(glycan_list[glycan_index]) == 10:
                    acet_cnt += 1
                    acet_flag.append(acet_cnt)
                    for k in range(0,5):
                        mol.gets_atoms(i,atom[k],prev_chend_id,acet_flag[i-1])
                    if i == units_num-1: #Changes to end monomers 
                        mol.atoms[-4][1] = "SN4"
                        mol.atoms[-4][4] = "S2"
                        mol.atoms[-4][7] = 54
                elif len(glycan_list[glycan_index]) == 12:
                    acet_cnt += 2
                    acet_flag.append(acet_cnt)
                    for k in range(0,6):
                        mol.gets_atoms(i,atom[k],prev_chend_id,acet_flag[i-1])
                    if i == units_num-1: #Changes to end monomers 
                        mol.atoms[-4][1] = "SN4"
                        mol.atoms[-4][4] = "S2"
                        mol.atoms[-4][7] = 54
                else:
                    print(glycan_list[glycan_index])
                    raise RuntimeError("More than 5 atoms in glycan for "\
                                       + str(glycan_index))
        prev_chend_id += bead_cnt

        if bondIn:
            prev_chend_id = 0
            for j in range(ncnf_bundles*ch_per_cnf):
                bead_cnt = 0
                for i in range(units_num):
                    glycan_index = i + (j-1)*units_num
                    bead_cnt += int(0.5*len(glycan_list[glycan_index])) 
                    if i == 0 and j == 0:
                        if len(glycan_list[glycan_index]) == 8:
                            for k in range(0,4):
                                mol.gets_bonds(i,bondIn[k],prev_chend_id,0,0)
                        elif len(glycan_list[glycan_index]) == 10:
                            for k in range(0,5):
                                mol.gets_bonds(i,bondIn[k],prev_chend_id,0,0)
                        elif len(glycan_list[glycan_index]) == 12:
                            for k in range(0,6):
                                mol.gets_bonds(i,bondIn[k],prev_chend_id,0,0)
                    else:
                        if len(glycan_list[glycan_index]) == 8:
                            for k in range(0,4):
                                mol.gets_bonds(i,bondIn[k],prev_chend_id,acet_flag[i-1],acet_flag[i-1])
                        elif len(glycan_list[glycan_index]) == 10:
                            for k in range(0,5):
                                mol.gets_bonds(i,bondIn[k],prev_chend_id,acet_flag[i-1],acet_flag[i-1]) 
                        elif len(glycan_list[glycan_index]) == 12:
                            for k in range(0,6):
                                mol.gets_bonds(i,bondIn[k],prev_chend_id,acet_flag[i-1],acet_flag[i-1]) 
                prev_chend_id += bead_cnt

        if bondConnect:
            prev_chend_id = 0
            for j in range(ncnf_bundles*ch_per_cnf):
                bead_cnt = 0
                for row in bondConnect:
                    for i in range(units_num - 1):
                        if i == 0 and j == 0:
                            mol.gets_bonds(i, row, prev_chend_id, 0, acet_flag[i])
                        else:
                            mol.gets_bonds(i, row, prev_chend_id, acet_flag[i-1],acet_flag[i])
                prev_chend_id += bead_cnt

        if angleIn:
            prev_chend_id = 0
            for j in range(ncnf_bundles*ch_per_cnf):
                bead_cnt = 0
                for i in range(units_num):
                    glycan_index = i + (j-1)*units_num
                    bead_cnt += int(0.5*len(glycan_list[glycan_index])) 
                    if len(glycan_list[glycan_index]) == 8:
                        for k in range(0,2):
                            mol.gets_angles(i, angleIn[k], prev_chend_id,0,0,0)
                    elif len(glycan_list[glycan_index]) == 10:
                        for k in range(0,3):
                            mol.gets_angles(i, angleIn[k], prev_chend_id,0,0,0)
                    elif len(glycan_list[glycan_index]) == 12:
                        for k in range(0,4):
                            mol.gets_angles(i, angleIn[k], prev_chend_id,0,0,0)
                prev_chend_id += bead_cnt

        if angleConnect:
            for row in angleConnect:
                for i in range(units_num-1):
                    mol.gets_angles(i, row, prev_chend_id, 0, 0, 0) # CHECK

        if dihIn:
            prev_chend_id = 0
            for j in range(ncnf_bundles*ch_per_cnf):
                bead_cnt = 0
                for i in range(units_num):
                    glycan_index = i + (j-1)*units_num
                    bead_cnt += int(0.5*len(glycan_list[glycan_index])) 
                    if len(glycan_list[glycan_index]) == 8:
                        for k in range(len(dihIn)-1):
                            mol.gets_dihedrals(i, dihIn[k], prev_chend_id,0,0,0,0)
                    elif len(glycan_list[glycan_index]) == 10:
                        for k in range(len(dihIn)):
                            mol.gets_dihedrals(i, dihIn[k], prev_chend_id,0,0,0,0)
                prev_chend_id += bead_cnt

        if dihConnect:
            for row in dihConnect:
                for i in range(units_num-1):
                    mol.gets_dihedrals(i, row, prev_chend_id, 0, 0, 0, 0) #CHECK

        if improperIn:
            prev_chend_id = 0
            for j in range(ncnf_bundles*ch_per_cnf):
                bead_cnt = 0
                for i in range(units_num):
                    glycan_index = i + (j-1)*units_num
                    bead_cnt += int(0.5*len(glycan_list[glycan_index])) 
                    if len(glycan_list[glycan_index]) == 8:
                        for k in range(len(improperIn)-1):
                            mol.gets_impropers(i, improperIn[k], prev_chend_id,0,0,0,0)
                    elif len(glycan_list[glycan_index]) == 10:
                        for k in range(len(improperIn)):
                            mol.gets_impropers(i, improperIn[k], prev_chend_id,0,0,0,0)
                prev_chend_id += bead_cnt

        if improperConnect:
            for row in improperConnect:
                for i in range(units_num-1):
                    mol.gets_impropers(i, row, prev_chend_id, 0, 0, 0, 0) #CHECK

    
    mol.writeTop(outfname,molname)
#----------------------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#----------------------------------------------------------------------------------
