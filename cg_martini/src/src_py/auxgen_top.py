# Auxiliary file for gen_top.py

import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
import struct
import pdb
#----------------------------------------------------------------
# Generate log file
def gen_logfile(fname,ncnf,acetfrac,cell_dp,ch_per_cnf):
    fout = open('log.dat','w')
    fout.write('Input file name: %s\n' %(fname))
    fout.write('Number of CNF chains: %d\n' %(ncnf))
    fout.write('Acetylated fraction: %g\n' %(acetfrac))
    fout.write('Cellulose degree of polymerization: %d\n' %(cell_dp))
    fout.write('Number of CNF chains per fiber: %d\n' %(ch_per_cnf))
    return fout
#----------------------------------------------------------------
# Read formatted GRO file
# Ref: https://github.com/hernanchavezthielemann/GRO2LAM/blob/27ene19/lib/handling/gromacs.py
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
            elif len(line.strip().split()) == 3:
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
#----------------------------------------------------------------
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
#----------------------------------------------------------------
# Create bond list between the following glycan type names
# O2-C4; C4-C6; C4-C4
def create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list):
    # Create empty list
    bond_list = [[] for i in range((3*cell_dp-1)*ncnf*ch_per_cnf)]
    bnd_cntr = 0
    fbnd = open('bond.txt','w')
    # Run through glycan list and connect "inner" beads (O2-C4;C4-C6)
    for i in range(len(glycan_list)):
        for j in range(0,len(glycan_list[i]),2):
            for k in range(j+2,len(glycan_list[i]),2):
                if (glycan_list[i][j+1] == 'O2' and glycan_list[i][k+1] == 'C4')\
                   or (glycan_list[i][j+1] == 'C4' and glycan_list[i][k+1] == 'C6'):
                    bond_list[bnd_cntr].append(glycan_list[i][j])
                    bond_list[bnd_cntr].append(glycan_list[i][k])
                    bnd_cntr += 1
                    fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
    # Connection between glycans (C4-C4)
    for ch in range(ncnf*ch_per_cnf):
        for glyres in range(cell_dp-1):
            i = cell_dp*ch + glyres
            j = glycan_list[i].index('C4') - 1; # -1 to point to ID of C4
            k = glycan_list[i+1].index('C4') - 1
            bond_list[bnd_cntr].append(glycan_list[i][j])
            bond_list[bnd_cntr].append(glycan_list[i+1][k])
            bnd_cntr += 1
            fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i+1][k]))
    return bond_list
#----------------------------------------------------------------
