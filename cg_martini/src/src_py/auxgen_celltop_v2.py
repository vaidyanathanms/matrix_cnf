# To generate itp files for certain polymers and acetylated cellulose
# Uses Martini V2.0 for cellulose structure
# Main file: gen_top.py
# Reference: Martini3-Polyply
# https://github.com/ricalessandri/Martini3-small-molecules/blob/main/models/martini_v3.0.0_small_molecules_v1.itp
# NOTE: Angle constants for PLA/P3HB end groups are assigned same value as those in the middle. Requires more work
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
#----------------------------------------------------------------------------------
# Bonded parameters between P6, PX and P5 are from Lopez et al. 2015 Cellulose paper based on martini v2.0
# P6 below (based on martini v3.0) is same as P4 in martini v2.0
# P5 below (based on martini v3.0) is same as P1 in martini v2.0 
# Bonded parameters between P5 and N4a are taken from lipids file
#----------------------------------------------------------------------------------
# Create bond list between the following glycan type names
def create_bond_list(cell_dp,ncnf_per_bundle,ch_per_cnf,glycan_list):
    # Create empty list
    bond_list = [[]]
    bnd_cntr = 0
    fbnd = open('bond.txt','w')
    # Run through glycan list and connect "inner" beads (P6-PX;PX-P5)
    for i in range(len(glycan_list)):
        for j in range(0,len(glycan_list[i]),2):
            for k in range(j+2,len(glycan_list[i]),2):
                if (glycan_list[i][j+1] == 'P6' and glycan_list[i][k+1] == 'PX'):  
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bond_list[bnd_cntr].append(1)
                   bond_list[bnd_cntr].append(0.23)
                   bond_list[bnd_cntr].append(30000)
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])
                elif (glycan_list[i][j+1] == 'PX' and glycan_list[i][k+1] == 'P5'):  
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bond_list[bnd_cntr].append(1)
                   bond_list[bnd_cntr].append(0.22)
                   bond_list[bnd_cntr].append(30000)
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])
                elif (glycan_list[i][j+1] == 'P5' and glycan_list[i][k+1] == 'SN4a'): # This is based on P1-Na bond parameters in martini_v2.0_lipids_all_201506.itp file
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bond_list[bnd_cntr].append(1)
                   bond_list[bnd_cntr].append(0.40)
                   bond_list[bnd_cntr].append(30000)
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])
                elif (glycan_list[i][j+1] == 'SN4a' and glycan_list[i][k+1] == 'TC1'):
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bond_list[bnd_cntr].append(1)
                   bond_list[bnd_cntr].append(0.37)
                   bond_list[bnd_cntr].append(1250)
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])

    # Connection between glycans (PX-PX)
    for ch in range(ncnf_per_bundle*ch_per_cnf):
        for glyres in range(cell_dp-1):
            i = cell_dp*ch + glyres
            j = glycan_list[i].index('PX') - 1; # -1 to point to ID of PX
            k = glycan_list[i+1].index('PX') - 1
            bond_list[bnd_cntr].append(glycan_list[i][j])
            bond_list[bnd_cntr].append(glycan_list[i+1][k])
            bond_list[bnd_cntr].append(1)
            bond_list[bnd_cntr].append(0.55)
            bond_list[bnd_cntr].append(30000)
            bnd_cntr += 1
            fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i+1][k]))
            bond_list.append([])
    del(bond_list[-1])
    return bond_list
#----------------------------------------------------------------------------------
# Create angle list between the following glycan type names
# P6-PX-P5; P6-PX-PX; P5-PX-PX; PX-PX-PX; PX-P5-SN4a, P5-SN4a-TC1  
def create_angle_list(cell_dp,ncnf_per_bundle,ch_per_cnf,glycan_list):
     # Create empty list
     angle_list = [[]] 
     angle_cntr = 0
     fang = open('angle.txt','w')
     # Run through glycan list and make angle connections - P6-PX-P5
     for i in range(len(glycan_list)):
         for j in range(0,len(glycan_list[i]),2):
             for k in range(j+2,len(glycan_list[i]),2):
                 for l in range(k+2,len(glycan_list[i]),2):
                     if (glycan_list[i][j+1] == 'P6' and glycan_list[i][k+1] == 'PX' and glycan_list[i][l+1] == 'P5'):
                        angle_list[angle_cntr].append(glycan_list[i][j])
                        angle_list[angle_cntr].append(glycan_list[i][k])
                        angle_list[angle_cntr].append(glycan_list[i][l])
                        angle_list[angle_cntr].append(2)
                        angle_list[angle_cntr].append(166.0)
                        angle_list[angle_cntr].append(50.0)
                        angle_cntr += 1
                        fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i][l]))
                        angle_list.append([])
     # Run through glycan list and make angle connections - PX-P5-SN4a
     for i in range(len(glycan_list)):
         for j in range(0,len(glycan_list[i]),2):
             for k in range(j+2,len(glycan_list[i]),2):
                 for l in range(k+2,len(glycan_list[i]),2):
                     if (glycan_list[i][j+1] == 'PX' and glycan_list[i][k+1] == 'P5' and glycan_list[i][l+1] == 'SN4a'):
                        angle_list[angle_cntr].append(glycan_list[i][j])
                        angle_list[angle_cntr].append(glycan_list[i][k])
                        angle_list[angle_cntr].append(glycan_list[i][l])
                        angle_list[angle_cntr].append(2)
                        angle_list[angle_cntr].append(100.0) # taken from DPMG lipid for P1-P1-Na martini_v2.0_lipids_all_201506.itp file 
                        angle_list[angle_cntr].append(35.0) # taken from DPMG lipid for P1-P1-Na martini_v2.0_lipids_all_201506.itp file
                        angle_cntr += 1
                        fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i][l]))
                        angle_list.append([])
     # Run through glycan list and make angle connections - P5-SN4a-TC1
     for i in range(len(glycan_list)):
         for j in range(0,len(glycan_list[i]),2):
             for k in range(j+2,len(glycan_list[i]),2):
                 for l in range(k+2,len(glycan_list[i]),2):
                     if (glycan_list[i][j+1] == 'P5' and glycan_list[i][k+1] == 'SN4a' and glycan_list[i][l+1] == 'TC1'):
                        angle_list[angle_cntr].append(glycan_list[i][j])
                        angle_list[angle_cntr].append(glycan_list[i][k])
                        angle_list[angle_cntr].append(glycan_list[i][l])
                        angle_list[angle_cntr].append(2)
                        angle_list[angle_cntr].append(180.0) # taken from LPC (IPC) lipid for P1-Na-C1 martini_v2.0_lipids_all_201506.itp file
                        angle_list[angle_cntr].append(25.0) # taken from LPC (IPC) lipid for P1-Na-C1 martini_v2.0_lipids_all_201506.itp file
                        angle_cntr += 1
                        fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i][l]))
                        angle_list.append([])

    # Run through glycan list and make angle connections - P6-PX-PX; P5-PX-PX for residues in the direction or order B1-B2-B5, B3-B2-B5            
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('P6') - 1; # -1 to point to ID of P6
             k = glycan_list[i].index('PX') - 1
             l = glycan_list[i+1].index('PX') - 1
             angle_list[angle_cntr].append(glycan_list[i][j])
             angle_list[angle_cntr].append(glycan_list[i][k])
             angle_list[angle_cntr].append(glycan_list[i+1][l])
             angle_list[angle_cntr].append(2)
             angle_list[angle_cntr].append(85.0)
             angle_list[angle_cntr].append(50.0)
             angle_cntr += 1
             fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i+1][l]))
             angle_list.append([])
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('P5') - 1; # -1 to point to ID of P5
             k = glycan_list[i].index('PX') - 1
             l = glycan_list[i+1].index('PX') - 1
             angle_list[angle_cntr].append(glycan_list[i][j])
             angle_list[angle_cntr].append(glycan_list[i][k])
             angle_list[angle_cntr].append(glycan_list[i+1][l])
             angle_list[angle_cntr].append(2)
             angle_list[angle_cntr].append(85.0)
             angle_list[angle_cntr].append(50.0)
             angle_cntr += 1
             fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i+1][l]))
             angle_list.append([])
     # Run through glycan list and make angle connections - P6-PX-PX; P5-PX-PX for residues in the direction B4-B5-B2, B6-B5-B2
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('P6') - 1; # -1 to point to ID of P6
             k = glycan_list[i].index('PX') - 1
             if i >= 0:
                l = glycan_list[i-1].index('PX') - 1
                angle_list[angle_cntr].append(glycan_list[i][j])
                angle_list[angle_cntr].append(glycan_list[i][k])
                angle_list[angle_cntr].append(glycan_list[i-1][l])
                angle_list[angle_cntr].append(2)
                angle_list[angle_cntr].append(80.0)
                angle_list[angle_cntr].append(50.0)

                angle_cntr += 1
                fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i-1][l]))
                angle_list.append([])
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('P5') - 1; # -1 to point to ID of P5
             k = glycan_list[i].index('PX') - 1
             if i >= 0:
                l = glycan_list[i-1].index('PX') - 1
                angle_list[angle_cntr].append(glycan_list[i][j])
                angle_list[angle_cntr].append(glycan_list[i][k])
                angle_list[angle_cntr].append(glycan_list[i-1][l])
                angle_list[angle_cntr].append(2)
                angle_list[angle_cntr].append(113.0)
                angle_list[angle_cntr].append(80.0)
                angle_cntr += 1
                fang.write('%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i-1][l]))
                angle_list.append([])
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('PX') - 1;
             k = glycan_list[i-1].index('PX') - 1
             if glyres > 0:
                l = glycan_list[i+1].index('PX') - 1
                angle_list[angle_cntr].append(glycan_list[i-1][k])
                angle_list[angle_cntr].append(glycan_list[i][j])
                angle_list[angle_cntr].append(glycan_list[i+1][l])
                angle_list[angle_cntr].append(2)
                angle_list[angle_cntr].append(165.0)
                angle_list[angle_cntr].append(50.0)
                angle_cntr += 1
                fang.write('%d\t%d\t%d\n' %(glycan_list[i-1][k],glycan_list[i][j],glycan_list[i+1][l]))
                angle_list.append([])
     del(angle_list[-1])
     return angle_list
#----------------------------------------------------------------------------------
# Create dihedral list between the following glycan type names
# P6-PX-PX-P6; P5-PX-PX-P5
def create_dihedral_list(cell_dp,ncnf_per_bundle,ch_per_cnf,glycan_list):
     # Create empty list
     dihed_list = [[]] 
     dihed_cntr = 0
     fdihed = open('dihedral.txt','w')
     # Run through glycan list and make angle connections - P6-PX-P5
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('P6') - 1;
             k = glycan_list[i].index('PX') - 1
             l = glycan_list[i+1].index('PX') - 1
             m = glycan_list[i+1].index('P6') - 1
             dihed_list[dihed_cntr].append(glycan_list[i][j])
             dihed_list[dihed_cntr].append(glycan_list[i][k])
             dihed_list[dihed_cntr].append(glycan_list[i+1][l])
             dihed_list[dihed_cntr].append(glycan_list[i+1][m])
             dihed_list[dihed_cntr].append(1)
             dihed_list[dihed_cntr].append(0.00)                
             dihed_list[dihed_cntr].append(10.00)
             dihed_list[dihed_cntr].append(1)
             dihed_cntr += 1
             fdihed.write('%d\t%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i+1][l],glycan_list[i+1][m]))
             dihed_list.append([])
     for ch in range(ncnf_per_bundle*ch_per_cnf):
         for glyres in range(cell_dp-1):
             i = cell_dp*ch + glyres
             j = glycan_list[i].index('P5') - 1;
             k = glycan_list[i].index('PX') - 1
             l = glycan_list[i+1].index('PX') - 1
             m = glycan_list[i+1].index('P5') - 1
             dihed_list[dihed_cntr].append(glycan_list[i][j])
             dihed_list[dihed_cntr].append(glycan_list[i][k])
             dihed_list[dihed_cntr].append(glycan_list[i+1][l])
             dihed_list[dihed_cntr].append(glycan_list[i+1][m])
             dihed_list[dihed_cntr].append(1)
             dihed_list[dihed_cntr].append(0.00) 
             dihed_list[dihed_cntr].append(10.00)
             dihed_list[dihed_cntr].append(1)
             dihed_cntr += 1
             fdihed.write('%d\t%d\t%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k],glycan_list[i+1][l],glycan_list[i+1][m]))
             dihed_list.append([])
     del(dihed_list[-1])
     return dihed_list
#----------------------------------------------------------------------------------------
def write_celltop(posre_fname,residarr,resnamearr,aidarr,anamearr,massarr,bond_list,angle_list,dihed_list):
     ftop = open('CELLULOSE_acetylated.top','w') 
     charge = 0.0
     rem_str = "; qtot 0"
     ftop.write(";\tFile CELLULOSE_acetylated.top was generated from gen_top.py\n")
     ftop.write(";\tInclude forcefield parameters\n")
#     ftop.write('#include "./all_toppar/forcefield.itp"\n\n') 
     ftop.write("[ moleculetype ]\n")
     ftop.write("; Name             nrexcl\n")
     ftop.write("CELLULOSE-acetylated       1\n\n") # only 1-2 bonded neighbors are excluded from LJ calculations in MARTINI
     ftop.write("[ atoms ]\n")
     ftop.write(";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n")
     for i in range(0, len(aidarr)):
         ftop.write("%d\t%s\t%d\t%s\t%s\t%d\t%f\t%s\n" %(i+1,anamearr[i],residarr[i],resnamearr[i],anamearr[i],residarr[i],charge,rem_str))           
     ftop.write("\n[ bonds ]\n") 
     for i in range(0, len(bond_list)): 
         ftop.write("%d\t%d\t%d\t%f\t%f\n" %(bond_list[i][0],bond_list[i][1],bond_list[i][2],bond_list[i][3],bond_list[i][4]))
     ftop.write("\n[ angles ]\n")
     for i in range(0, len(angle_list)):
         ftop.write("%d\t%d\t%d\t%d\t%f\t%f\n" %(angle_list[i][0],angle_list[i][1],angle_list[i][2],angle_list[i][3],angle_list[i][4],angle_list[i][5]))
     ftop.write("\n[ dihedrals ]\n")
     for i in range(0, len(dihed_list)):
         ftop.write("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\n\n" %(dihed_list[i][0],dihed_list[i][1],dihed_list[i][2],dihed_list[i][3],dihed_list[i][4],dihed_list[i][5],dihed_list[i][6],dihed_list[i][7]))
     ftop.write("; Include Position restraint file\n")
     ftop.write("#ifdef POSRES\n")
     ftop.write("#include \"posre-%s.itp\"\n" %(posre_fname))
     ftop.write("#endif\n")
     return ftop   
#----------------------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#----------------------------------------------------------------------------------
