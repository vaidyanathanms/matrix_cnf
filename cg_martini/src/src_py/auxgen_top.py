# To generate itp files for certain polymers and acetylated cellulose
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
# Make output directory
def make_out_dir(currdir,moltype):
    pathstrs = currdir.split("/")
    parpath  = "/".join( str(x) for x in pathstrs[:len(pathstrs)-2])
    outpath  = parpath+'/'+moltype
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    return outpath
#----------------------------------------------------------------------------------
# Generate log file for cellulose/acetylated cellulose systems
def gen_logfile(fname,ncnf,acetfrac,cell_dp,ch_per_cnf):
    fout = open('log.dat','w')
    fout.write('Input file name: %s\n' %(fname))
    fout.write('Number of CNF chains: %d\n' %(ncnf))
    fout.write('Acetylated fraction: %g\n' %(acetfrac))
    fout.write('Cellulose degree of polymerization: %d\n' %(cell_dp))
    fout.write('Number of CNF chains per fiber: %d\n' %(ch_per_cnf))
    return fout
#-----------------------------------------------------------------------------------
# Read formatted GRO file
# Ref: https://github.com/hernanchavezthielemann/GRP4LAM/blob/27ene19/lib/handling/gromacs.py
def read_gro_file(inpfyle):
    fmt_at = b'5s5s5s5s8s8s8s'
    aidarr = []; anamearr = [];residarr = []; resnamearr = []
    rxarr = []; ryarr = []; rzarr = []; massarr = []
    with open(inpfyle,'r') as fin:
        fin.readline()
        nbeads = int(fin.readline().strip())
        for line in fin:
            if line.startswith(';'): # skip lines starting with #
                continue
            if not line: # skip empty lines
                continue
            elif len(line.strip().split()) == 3 or len(line.strip().split()) == 9:
                lbox = line.strip().split()
            else:
                residarr.append(int(line[:5].strip()))
                resnamearr.append(line[5:10].strip())
                anamearr.append(line[10:15].strip())
                aidarr.append(int(line[15:20].strip()))
                rxarr.append(float(line[20:28].strip()))
                ryarr.append(float(line[28:36].strip()))
                rzarr.append(float(line[36:44].strip()))

            if 'H' in line[10:15].strip():
                mval = 1.008
            elif 'C' in line[10:15].strip():
                mval = 12.0108
            elif 'O' in line[10:15].strip():
                mval = 15.9994
            elif 'S' in line[10:15].strip():
                mval = 32.065
            elif 'N' in line[10:15].strip():
                mval = 14.0067
            else: # Default Martini
                mval = 72
            massarr.append(mval)
    if len(residarr) == 0:
        raise RuntimeError('No atoms read - file corrupted?')
    return residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr
#----------------------------------------------------------------------------------
# Creating MARTINI beads list
def create_martini_beads(cell_dp,ncnf,ch_per_cnf,residarr,\
                         aidarr,anamearr):
    # Create empty array
    glycan_list = [[] for i in range(cell_dp*ncnf*ch_per_cnf)]
    # Sort beads
    headid_ptr = -1; resid_gro = 0; glycan_cnt = -1
    while resid_gro < len(residarr):
        if residarr[resid_gro] == headid_ptr: #check if resID is same
            while residarr[resid_gro] == headid_ptr:
                glycan_list[glycan_cnt].append(aidarr[resid_gro])
                glycan_list[glycan_cnt].append(anamearr[resid_gro])
                resid_gro += 1 #move to next "atom" in gro file
                if resid_gro == len(residarr): #break while loop
                    break
        else: #update new glycan
            glycan_cnt += 1 #move to new "glycan"
            headid_ptr = residarr[resid_gro] #update pointer to new glycan
    return glycan_list
#----------------------------------------------------------------------------------
# Bonded parameters between P6, PX and P5 are from Lopez et al. 2015 Cellulose paper based on martini v2.0
# P6 below (based on martini v3.0) is same as P4 in martini v2.0
# P5 below (based on martini v3.0) is same as P1 in martini v2.0 
# Bonded parameters between P5 and N4a are taken from lipids file
#----------------------------------------------------------------------------------
# Create bond list between the following glycan type names
def create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list):
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
    for ch in range(ncnf*ch_per_cnf):
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
def create_angle_list(cell_dp,ncnf,ch_per_cnf,glycan_list):
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
     for ch in range(ncnf*ch_per_cnf):
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
     for ch in range(ncnf*ch_per_cnf):
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
     for ch in range(ncnf*ch_per_cnf):
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
     for ch in range(ncnf*ch_per_cnf):
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
     for ch in range(ncnf*ch_per_cnf):
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
def create_dihedral_list(cell_dp,ncnf,ch_per_cnf,glycan_list):
     # Create empty list
     dihed_list = [[]] 
     dihed_cntr = 0
     fdihed = open('dihedral.txt','w')
     # Run through glycan list and make angle connections - P6-PX-P5
     for ch in range(ncnf*ch_per_cnf):
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
     for ch in range(ncnf*ch_per_cnf):
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
def write_celltop(residarr,resnamearr,aidarr,anamearr,massarr,bond_list,angle_list,dihed_list):
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
         ftop.write("%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\n" %(dihed_list[i][0],dihed_list[i][1],dihed_list[i][2],dihed_list[i][3],dihed_list[i][4],dihed_list[i][5],dihed_list[i][6],dihed_list[i][7]))
     return ftop   
#----------------------------------------------------------------------------------
# Design PETG
def design_petg(nmons,moltype):

    # Ref: See methyl-benzoate in the main reference (Martini3-Polyply) on headers
    # Structure
    # TP2-TC5-TC5-N4a-N4a-TC5-TC5-N4a-...TC5-TC5-N4a-TP1
    #   \  \ /  /       \  \ /  /         \ /   /
    #      TC5             TC5            TC5

    natoms = 1 + 5*nmons; nbonds = 5*nmons; nconst = 3*nmons;
    ndihds = 2*nmons; nexcls = 2*nmons

    at_list = [[0 for x in range(4)] for y in range(natoms)] #atoms
    bo_list = [[0 for x in range(5)] for y in range(nbonds)] #bonds
    co_list = [[0 for x in range(4)] for y in range(nconst)] #constraints
    di_list = [[0 for x in range(7)] for y in range(ndihds)] #dihedrals
    ex_list = [[0 for x in range(2)] for y in range(nexcls)] #exclusions
    
    for atcnt in range(natoms):
        at_list[atcnt][0] = atcnt+1
        
    atcnt = 0
    for mon_index in range(1,nmons+1):

        if mon_index == 1:

            for j in range(5):
                if j == 0: 
                    at_list[atcnt+j][1] = 'TP2' #Terminal acid
                    at_list[atcnt+j][2] = 'TP2'
                    at_list[atcnt+j][3] = 36.0
                elif j == 4:
                    at_list[atcnt+j][1] = 'N4a' #Ester
                    at_list[atcnt+j][2] = 'N4a'
                    at_list[atcnt+j][3] = 72.0
                else:
                    at_list[atcnt+j][1] = 'TC5' # Ring
                    at_list[atcnt+j][2] = 'TC5'
                    at_list[atcnt+j][3] = 36.0


            atcnt = atcnt + 5

        else:
            for j in range(5):
                if j == 0 or j == 4:
                    at_list[atcnt+j][1] = 'N4a'
                    at_list[atcnt+j][2] = 'N4a'
                    at_list[atcnt+j][3] = 72.0
                else:
                    at_list[atcnt+j][1] = 'TC5'
                    at_list[atcnt+j][2] = 'TC5'
                    at_list[atcnt+j][3] = 36.0
                    
            atcnt = atcnt + 5

    at_list[atcnt][1] = 'TP1' 
    at_list[atcnt][2] = 'TP1' #Terminal alcohol
    at_list[atcnt][3] = 36.0
    atcnt += 1 # last atom
    
    # Bond list
    bcnt = 0
    for mon_index in range(1,nmons+1):
        fatom = 1+5*(mon_index-1) #first atom in the monomer
        for j in range(5):
            bo_list[bcnt][0] = fatom+j
            bo_list[bcnt][1] = fatom+j+1
            bo_list[bcnt][2] = 1 # Bond type - change if needed
            if j == 0: #"Left end" of a monomer (terminal TP2 or N4a)
                #Bond with adjacent TC5 
                bo_list[bcnt][3] = 0.389 #bond length
                bo_list[bcnt][4] = 5000 #k_b
                #Bond with other TC5
                bo_list[bcnt+1][0] = fatom+j
                bo_list[bcnt+1][1] = fatom+j+2
                bo_list[bcnt+1][2] = 1
                bo_list[bcnt+1][3] = 0.389 #bond length
                bo_list[bcnt+1][4] = 5000 #k_b
                bcnt += 2
            elif j == 3: #"Right end" of a ring (TC5)
                # Bond with next N4a
                bo_list[bcnt][3] = 0.389 #bond length
                bo_list[bcnt][4] = 5000 #k_b
                bcnt += 1
            elif j == 4: #"Right end" of a monomer (N4a)
                # Bond with next N4a/terminal TP1
                bo_list[bcnt][3] = 0.37 #bond length
                bo_list[bcnt][4] = 1250 #k_b
                #Bond with previous TC5
                bo_list[bcnt+1][0] = fatom+j
                bo_list[bcnt+1][1] = fatom+j-2
                bo_list[bcnt+1][2] = 1
                bo_list[bcnt+1][3] = 0.37 #bond length
                bo_list[bcnt+1][4] = 1250 #k_b
                bcnt += 2

    # Constraint list
    ccnt = 0
    for mon_index in range(1,nmons+1):
        f_ringid = 5*(mon_index-1)+2 #first ring atom id
        for j in range(2):
            co_list[ccnt][0] = f_ringid+j
            co_list[ccnt][1] = f_ringid+j+1 
            co_list[ccnt][2] = 1 #function
            co_list[ccnt][3] = 0.264 #constraint length
            ccnt += 1
        # Extra bond between "first" and "last" ring atom
        co_list[ccnt][0] = f_ringid
        co_list[ccnt][1] = f_ringid+2 
        co_list[ccnt][2] = 1 #function
        co_list[ccnt][3] = 0.264 #constraint length
        ccnt += 1

    # Dihedral/exclusion list
    dicnt = 0; excnt = 0
    for mon_index in range(1,nmons+1):
        fatom = 1 + 5*(mon_index-1) #first atom in the monomer
        for j in range(2):
            di_list[dicnt][0] = fatom + j
            di_list[dicnt][1] = fatom + j + 1
            di_list[dicnt][2] = fatom + j + 2 
            di_list[dicnt][3] = fatom + j + 3
            di_list[dicnt][4] = 2 # function
            di_list[dicnt][5] = 180.0 #angle
            di_list[dicnt][6] = 50.0 #force constant
            ex_list[excnt][0] = fatom + j
            ex_list[excnt][1] = fatom + j + 3
            dicnt += 1; excnt += 1

    # Sanity check
    sanity_check(natoms,atcnt,nbonds,bcnt,0,0,ndihds,dicnt,nconst,ccnt)

    # Write headers
    outdir   = make_out_dir(os.getcwd(),moltype.lower())
    fpetg = open(outdir+'/'+moltype.lower()+'_'+str(nmons)+'.itp','w')
    write_header(fpetg,nmons,natoms,nbonds,0,ndihds,nconst,moltype) #headers
    fpetg.write('; Ref1: Martini V2.0 Lipids for N4a-N4a\n')
    fpetg.write('; Ref2: Martini V3.0 for others\n')
    
    # Write topology
    write_topology(fpetg,moltype,nmons,natoms,nbonds,nconst,0,ndihds,nexcls,\
                   at_list,bo_list,co_list,[],di_list,ex_list)

    # Close file
    fpetg.close()
#----------------------------------------------------------------------------------
# Design PLA
def design_pla(nmons,moltype):

    # Ref: MartiniV2.0-Lipids and Li et al., Acta mechanica solida sinica 30 (2017) 630–637
    # Structure
    # TP2-N4a-N4a-N4a-N4a-....-N4a-TC1
    # The structure in atomistic is C(OH)O-(CH(CH3)-C(O))n-CH2-CH3
    # NOTE: The monomer of PLA is different from what is defined as
    # N4a. Number of N4a CG beads = nmons-1
    # Also, note that the size of the monomer is changed to
    # incorporate the extra CH3 (see figure below).
    # N4a -> -C-O-CH-
    #         ||   |
    #         O   CH3
    natoms = 1 + nmons; nbonds = natoms-1; nangls = natoms-2
    ndihds = natoms-3; nexcls = 2*nmons

    at_list = [[0 for x in range(4)] for y in range(natoms)] #atoms
    bo_list = [[0 for x in range(5)] for y in range(nbonds)] #bonds
    an_list = [[0 for x in range(6)] for y in range(nbonds)] #angles
    di_list = [[0 for x in range(7)] for y in range(ndihds)] #dihedrals
    ex_list = [[0 for x in range(2)] for y in range(nexcls)] #exclusions
    
    for atcnt in range(natoms):
        at_list[atcnt][0] = atcnt+1

    # Atom list
    for at_index in range(natoms):
        if at_index == 0:
            at_list[at_index][1] = 'TP2'  # COOH Terminal end
            at_list[at_index][2] = 'TP2'
            at_list[at_index][3] = 36.0
        elif at_index == natoms-1:
            at_list[at_index][1] = 'TC1' # CH2-CH3 Terminal end
            at_list[at_index][2] = 'TC1'
            at_list[at_index][3] = 36.0
        else:
            at_list[at_index][1] = 'N4a' # Modified N4a 
            at_list[at_index][2] = 'N4a'
            at_list[at_index][3] = 72.0

    # Bond list
    for bt_index in range(nbonds):
        bo_list[bt_index][0] = bt_index+1
        bo_list[bt_index][1] = bt_index+2
        bo_list[bt_index][2] = 1 # Bond type - change if needed
        if bt_index == 0: #Terminal TP2-N4a
            bo_list[bt_index][3] = 0.5 #bond length
            bo_list[bt_index][4] = 10000 #k_b
        elif bt_index == nbonds-1: #Terminal N4a-TC1
            bo_list[bt_index][3] = 0.47 #bond length
            bo_list[bt_index][4] = 3600 #k_
        else: #N4a-N4a (From Li et al., see above)
            bo_list[bt_index][3] = 0.434  #bond length or 0.37 in
            #Martini V2.0
            bo_list[bt_index][4] = 8000 #k_b or 1250 in Martini V2.0

    # Angle list
    for an_index in range(nangls):
        an_list[an_index][0] = an_index+1
        an_list[an_index][1] = an_index+2
        an_list[an_index][2] = an_index+3
        an_list[an_index][3] = 2 # Angle type - change if needed
        if an_index == 0: #Terminal TC1-N4a-N4a
            an_list[an_index][4] = 120 #eqbm angle
            an_list[an_index][5] = 118 #k_theta
        elif an_index == nangls-1: #Terminal N4a-N4a-TP2
            an_list[an_index][4] = 120 #eqbm angle
            an_list[an_index][5] = 118 #k_theta
        else: #N4a-N4a-N4a (Li et al., Acta mechanica solida sinica 30 (2017) 630–637)
            an_list[an_index][4] = 120 #eqbm angle
            an_list[an_index][5] = 118 #k_theta

    # Sanity check
    sanity_check(natoms,1+at_index,nbonds,1+bt_index,nangls,1+an_index,\
                 0,0,0,0)

    # Write headers
    outdir   = make_out_dir(os.getcwd(),moltype.lower())
    fpla = open(outdir+'/'+moltype.lower()+'_'+str(nmons)+'.itp','w')
    write_header(fpla,nmons,natoms,nbonds,nangls,0,0,moltype) #headers
    fpla.write('; Ref1: Martini V2.0 Lipids\n')
    fpla.write('; Ref2: for angles/bonds-Li et al., Acta mechanica solida '+\
               'sinica 30 (2017) 630–637\n')
    
    # Write topology
    write_topology(fpla,moltype,nmons,natoms,nbonds,0,nangls,0,0,\
                   at_list,bo_list,[],an_list,[],[])

    # Close file
    fpla.close()
#----------------------------------------------------------------------------------      
# Design P3HB
def design_p3hb(nmons,moltype):

    # Ref: Martini V2.0 Lipids and Rossi et al., Macromolecules, 2011,
    # 44, 15, 6198-6208
    # Structure
    # TP2-C2-TN4a-C2-TN4a-...-C2-TN4a-SC2
    # NOTE: The mapping reduces the definition of monomer by one. So
    # if there are n-monomers, there will be 2+2*(n-1) beads, where
    # each mapped monomer is defined by two beads - C2-TN4a
    # TN4a ->   -O-C-    C2 -> CH2-CH-
    #             ||                |
    #             O                CH3
    natoms = 2 + 2*(nmons-1); nbonds = natoms-1; nangls = natoms-2 
    ndihds = natoms-3; nexcls = 2*nmons

    at_list = [[0 for x in range(4)] for y in range(natoms)] #atoms
    bo_list = [[0 for x in range(5)] for y in range(nbonds)] #bonds
    an_list = [[0 for x in range(6)] for y in range(nangls)] #angles
    di_list = [[0 for x in range(7)] for y in range(ndihds)] #dihedrals
    ex_list = [[0 for x in range(2)] for y in range(nexcls)] #exclusions
    
    for atcnt in range(natoms):
        at_list[atcnt][0] = atcnt+1
        
    atcnt = 0
    # Atom list
    for at_index in range(natoms):
        if at_index == 0:
            at_list[at_index][1] = 'TP2' #Acid end
            at_list[at_index][2] = 'TP2'
            at_list[at_index][3] = 36.0
        elif at_index == natoms-1:
            at_list[at_index][1] = 'SC1' #alkane end
            at_list[at_index][2] = 'SC1'
            at_list[at_index][3] = 54.0
        else:
            if at_index%2 == 0: # Second atom in backbone sequence
                at_list[at_index][1] = 'TN4a' # O-C=O
                at_list[at_index][2] = 'TN4a' #Ester
                at_list[at_index][3] = 36.0
            else: # First atom in backbone sequence
                at_list[at_index][1] = 'SC2' # CH2-CH(CH3)-
                at_list[at_index][2] = 'SC2' # Alkyl
                at_list[at_index][3] = 72.0
                
    # Bond list
    for bt_index in range(nbonds):
        bo_list[bt_index][0] = bt_index+1
        bo_list[bt_index][1] = bt_index+2
        bo_list[bt_index][2] = 1 # Bond type - change if needed
        if bt_index == 0: #Terminal TP2-C2 (based on P2-C1 and C1-C1)
            bo_list[bt_index][3] = 0.47 #bond length
            bo_list[bt_index][4] = 1250 #k_b
        elif bt_index == nbonds-1: #Terminal TN4a-SC1 (Martini V2.0
            #for Na-C1)
            bo_list[bt_index][3] = 0.47 #bond length
            bo_list[bt_index][4] = 3600 #k_b
        else: #C2-TN4a (based on Martini V2.0)
            bo_list[bt_index][3] = 0.47 #bond length
            bo_list[bt_index][4] = 1250 #k_b

    # Angle list
    for an_index in range(nangls):
        an_list[an_index][0] = an_index+1
        an_list[an_index][1] = an_index+2
        an_list[an_index][2] = an_index+3
        an_list[an_index][3] = 1 # Angle type - change if needed
        if an_index == 0: #Terminal TP2-C2-TN4a
            an_list[an_index][4] = 120 #eqbm angle
            an_list[an_index][5] = 25 #k_theta
        elif an_index == nangls-1: #Terminal C2-TN4a-P1
            an_list[an_index][4] = 120 #eqbm angle
            an_list[an_index][5] = 25 #k_theta
        elif an_index%2 == 0: #TN4a-C2-TN4a Rossi, Macromolecules 2011
            an_list[an_index][4] = 130 #eqbm angle
            an_list[an_index][5] = 45 #k_theta   
        else: #C2-TN4a-C2 Rossi, Macromolecules 2011
            an_list[an_index][4] = 125 #eqbm angle
            an_list[an_index][5] = 25 #k_theta

    # Sanity check
    sanity_check(natoms,1+at_index,nbonds,1+bt_index,nangls,1+an_index,\
                 0,0,0,0)

    # Write headers
    outdir   = make_out_dir(os.getcwd(),moltype.lower())
    fp3hb = open(outdir+'/'+moltype.lower()+'_'+str(nmons)+'.itp','w')
    write_header(fp3hb,nmons,natoms,nbonds,nangls,0,0,moltype) #headers
    fp3hb.write('; Ref1:  Martini V2.0 Lipids\n')
    fp3hb.write('; Ref2: For angles Rossi et al., Macromolecules,'+\
                '2011, 44, 15, 6198-6208\n')
    
    # Write topology
    write_topology(fp3hb,moltype,nmons,natoms,nbonds,0,nangls,0,0,\
                   at_list,bo_list,[],an_list,[],[])

    # Close file
    fp3hb.close()
#----------------------------------------------------------------------------------
# Write headers
def write_header(fin,nmons,natoms,nbonds,nangls,ndihds,nconsts,moltype):
    fin.write('; itp file generated from scratch\n')
    fin.write('; Number of atomistic monomers per chain: %d\n' %(nmons))
    fin.write('; Number of CG beads per chain: %d\n' %(natoms))
    fin.write('; Number of CG bonds per chain: %d\n' %(nbonds))
    fin.write('; Number of CG angles per chain: %d\n' %(nangls))
    fin.write('; Number of CG dihedrals per chain: %d\n' %(ndihds))
    fin.write('; Number of CG constraints per chain: %d\n' %(nconsts))
    fin.write('\n')
    fin.write('[ moleculetype ]\n')
    fin.write('; name    nrexcl\n')
    fin.write(str(moltype) + '  1\n')
    fin.write('\n')
#----------------------------------------------------------------------------------
# Sanity cross check
def sanity_check(natoms,atcnt,nbonds,bcnt,nangls,ancnt,ndihds,dicnt,nconst,ccnt):
    if natoms != atcnt:
        raise RuntimeError('#atoms do not match:'+str(natoms)+' '+str(atcnt))
    if nbonds != bcnt:
        raise RuntimeError('#bonds do not match:'+str(nbonds)+' '+str(bcnt))
    if nangls != ancnt:
        raise RuntimeError('#bonds do not match:'+str(nangls)+' '+str(ancnt))
    if ndihds != dicnt:
        raise RuntimeError('#dihedrals do not match:'+str(ndihds)+' '+str(dicnt))
    if nconst != ccnt:
        raise RuntimeError('#constraints do not match:'+str(nconst)+' '+str(ccnt))
#----------------------------------------------------------------------------------
# Write topology to file
def write_topology(fin,moltype,nmons,natoms,nbonds,nconsts,nangls,ndihds,nexcls,\
                   at_list,bo_list,co_list,an_list,di_list,ex_list):

    if natoms != 0:
        #Write atoms
        fin.write('[ atoms ]\n')
        fin.write('; atomID moltype atomtyp atomname atomID charge mass\n')
        for i in range(natoms):
            fin.write('%6d %5s %6d %5s %5s %g %g %g\n' \
                      %(at_list[i][0],at_list[i][1],at_list[i][0],\
                        moltype,at_list[i][2],at_list[i][0],0.0,at_list[i][3]))
    else:
        raise RuntimeError('Unphysical system - zero atoms found')
    
    if nbonds != 0:
        #Write bonds
        fin.write('\n')
        fin.write('[ bonds ]\n')
        fin.write('; i   j    funct   length   force.c.\n')
        for i in range(nbonds):
            fin.write('%5d %5d %d %g %g\n' \
                      %(bo_list[i][0],bo_list[i][1],bo_list[i][2],\
                        bo_list[i][3],bo_list[i][4]))

    if nconsts != 0:
        #Write constraints
        fin.write('\n')
        fin.write('[ constraints ]\n')
        fin.write('; i   j    funct   length\n')
        for i in range(nconsts):
            fin.write('%5d %5d %g %g\n' \
                      %(co_list[i][0],co_list[i][1],\
                        co_list[i][2],co_list[i][3]))

    if nangls != 0:
        #Write angles
        fin.write('\n')
        fin.write('[ angles ]\n')
        fin.write('; i   j   k   funct   ref.angle  force_k\n')
        for i in range(nangls):
            fin.write('%5d %5d %5d %d %g %g\n' \
                      %(an_list[i][0],an_list[i][1],an_list[i][2],\
                        an_list[i][3],an_list[i][4],an_list[i][5]))

    if ndihds != 0:
        #Write dihedrals
        fin.write('\n')
        fin.write('[ dihedrals ]\n')
        fin.write('; i   j   k   l  funct   ref.angle  force_k\n')
        for i in range(ndihds):
            fin.write('%5d %5d %5d %5d %d %g %g\n' \
                      %(di_list[i][0],di_list[i][1],di_list[i][2],\
                        di_list[i][3],di_list[i][4],di_list[i][5],\
                        di_list[i][6]))
    if nexcls != 0:
        #Write exclusions
        fin.write('\n')
        fin.write('[ exclusions ]\n')
        fin.write('; i   j\n')
        for i in range(nexcls):
            fin.write('%5d %5d\n' %(ex_list[i][0],ex_list[i][1]))
#----------------------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#----------------------------------------------------------------------------------
