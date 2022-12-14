# To generate martini beads for unmodified/acetylated cellulose
# This file creates glycan_list and all generic i/p, o/p functions
# Main file: gen_top.py
# Can be used with auxgen_celltop_v2/v3.py
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

            if 'S' == line[10]: # small
                mval = 54
            elif 'T' == line[10]: # tiny
                mval = 36
            else: # Default Martini regular
                mval = 72
            massarr.append(mval)
    if len(residarr) == 0:
        raise RuntimeError('No atoms read - file corrupted?')
    return residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr
#----------------------------------------------------------------------------------
# Creating MARTINI beads list
# Returns glycan_list of size [nresidues][8] 
# 8 in second dimension comes from atomID and atomname per bead and
# there are 4 beads per residue
def create_martini_beads(cell_dp,ncnf_per_bundle,ch_per_cnf,residarr,\
                         aidarr,anamearr):
    # Create empty array
    glycan_list = [[] for i in range(cell_dp*ncnf_per_bundle*ch_per_cnf)]
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
