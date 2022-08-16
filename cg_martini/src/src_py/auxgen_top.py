# Auxiliary file for gen_top.py

import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
#----------------------------------------------------------------
# Generate log file
def gen_logfile(fname,ncnf,acetfrac,cell_dp,ch_per_cnf):
    fout = open('log.dat','w')
    fout.write('Input file name: ', fname)
    fout.write('Number of CNF chains: ', ncnf)
    fout.write('Acetylated fraction: ', acetfrac)
    fout.write('Cellulose degree of polymerization: ', cell_dp)
    fout.write('Number of CNF chains per fiber: ',ch_per_cnf)
    return fout
#----------------------------------------------------------------
# Creating MARTINI beads list
def create_martini_beads():
    # Create empty array
    glycan_list = [[] for i in range(cell_dp*ncnf*ch_per_cnf)]
    # Sort beads
    ref_glycan_id = -1; glycan_id = 0
    while glycan_id <= len(residarr):
        if residarr[glycan_id] != ref_glycan_id:
            ref_glycan_id = residarr[glycan_id]
            glycan_list[ref_glycan_id].append(aidarr[glycan_id])
            glycan_list[ref_glycan_id].append(anamearr[glycan_id])
            glycan_id += 1
        else:
            while residarr[glycan_id] != ref_glycan_id:
                glycan_list[ref_glycan_id].append(aidarr[glycan_id])
                glycan_list[ref_glycan_id].append(anamearr[glycan_id])
                glycan_id +=1
#----------------------------------------------------------------
# Read GRO file
def read_gro_file(inpfyle):
    fmt_at = b'5s5s5s8s8s8s'
    aidarr = []; anamearr = [];residarr = []; resnamearr = []
    rxarr = []; ryarr = []; rzarr = []; massarr = []
    with open(inpfyle,'r') as fin:
        fin.readline()
        nbeads = int(fin.readline.strip())
        for line in fin:
            if line.startswith(';'): # skip lines starting with #
                continue
            if not line: # skip empty lines
                continue
            elif len(line.strip.split()) == 3:
                lbox = line.strip.split():
            else:
                (resid,resname,aname,aid,rx,ry,rz)
                = struct.unpack(fmt_at,line.strip().encode('utf-8'))
   
            residarr.append(int(resid.decode('utf-8')))
            resnamearr.append((resname.decode('utf-8'))
            aidarr.append(int(aid.decode('utf-8')))
            anamearr.append(aname.decode('utf-8'))
            rxarr.append(float(rx.decode('utf-8')))
            ryarr.append(float(ry.decode('utf-8')))
            rzarr.append(float(rz.decode('utf-8')))

            if elem.decode('utf-8').strip() == 'H':
                mval = 1.008
            elif elem.decode('utf-8').strip() == 'C':
                mval = 12.0108
            elif elem.decode('utf-8').strip() == 'O':
                mval = 15.9994
            elif elem.decode('utf-8').strip() == 'S':
                mval = 32.065
            elif elem.decode('utf-8').strip() == 'N':
                mval = 14.0067
            else: # Default Martini
                mval = 72
            massarr.append(mval)

    return residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr
#----------------------------------------------------------------
