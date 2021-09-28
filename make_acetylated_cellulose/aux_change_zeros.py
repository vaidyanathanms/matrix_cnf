# Auxiliary file to convert zeros in input pdb files
import os
import sys
import numpy as np
import re
import shutil
import glob
import math
import struct
import subprocess
import random
#------------------------------------------------------------------

def read_and_process_file(inpfyle,outfyle):
    fmt_at = b'4sx6sx4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s'
    fmt_atq= b'4sx6sx4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s2s'
    fmt_ht = b'6s5s4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s'
    fmt_htq= b'6s5s4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s2s'

    recarr = []; massarr = []
    aidarr = []; anamearr = []
    residarr = [];resnamearr = []
    rxarr = []; ryarr = []; rzarr = []
    occarr = [];betaarr = [];elemarr = [];segnamearr=[]

    with open(inpfyle,'r') as fin:
        fin.readline()
        for line in fin:
            if line[0:4] == 'ATOM' and len(line.strip()) == 78:
                (rec,aid,aname,altloc,resname,chid,resid,rescode,\
                 rx,ry,rz,occ,beta,segname,elem) \
                    = struct.unpack(fmt_at,line.strip().encode('utf-8'))
                
            elif line[0:4] == 'ATOM' and len(line.strip()) == 80:
                (rec,aid,aname,altloc,resname,chid,resid,rescode,\
                 rx,ry,rz,occ,beta,segname,elem,qval) \
                    = struct.unpack(fmt_atq,line.strip().encode('utf-8'))
                qarr.append(qval.decode('utf-8'))
            elif line[0:4] == 'HEATOM' and len(line.strip()) == 78:
                (rec,aid,aname,altloc,resname,chid,resid,rescode,\
                 rx,ry,rz,occ,beta,segname,elem) \
                    = struct.unpack(fmt_ht,line.strip().encode('utf-8'))
            elif line[0:4] == 'HEATOM' and len(line.strip()) == 80:
                (rec,aid,aname,altloc,resname,chid,resid,rescode,\
                 rx,ry,rz,occ,beta,segname,elem,qval) \
                    = struct.unpack(fmt_htq,line.strip().encode('utf-8'))
                qarr.append(qval.decode('utf-8'))
            elif line[0:3] == 'END' or line[0:3] == 'TER':
                continue
            else:
                raise RuntimeError("Unknown line format", line)
            
            #atom details (altloc is not saved)
            recarr.append(rec.decode('utf-8'))
            aidarr.append(int(aid.decode('utf-8')))
            anamearr.append(aname.decode('utf-8'))

            #residue details (chid/rescode is not saved)
            resnamearr.append(resname.decode('utf-8'))
            residarr.append(int(resid.decode('utf-8')))

            #position details
            rxarr.append(float(rx.decode('utf-8')))
            ryarr.append(float(ry.decode('utf-8')))
            rzarr.append(float(rz.decode('utf-8')))

            #extra/segment details (qval is not saved)
            occarr.append(occ.decode('utf-8'))
            betaarr.append(beta.decode('utf-8'))
            elemarr.append(elem.decode('utf-8'))
            segnamearr.append(segname.decode('utf=8'))
            
            if elem.decode('utf-8').strip() == 'H':
                mval = 1.008
            elif elem.decode('utf-8').strip() == 'C':
                mval = 12.0108
            elif elem.decode('utf-8').strip() == 'O':
                mval = 15.9994

            massarr.append(mval)

    process_and_write(recarr,aidarr,anamearr,resnamearr,residarr\
                      ,rxarr,ryarr,rzarr,occarr,betaarr,elemarr,\
                      segnamearr,outfyle)
#------------------------------------------------------------------
# Find output write type
# ref: https://stackoverflow.com/questions/17796017/how-do-i-output-a-pdb-file-using-python-script
def process_and_write(recarr,aidarr,anamearr,resnamearr,residarr,
                      rxarr,ryarr,rzarr,occarr,betaarr,elemarr,\
                      segnamearr,outfyle):

    rdist  = 0.5 # move maximum by this distance

    fout = open(outfyle,'w')
    fout.write('%s\n' %('REMARK Processed output file VMS'))

    #NOTES: ignores altloc,chidarr,segname,qval
    #NOTES: In stackoverflow j[4] is defined which is ignored
    #j[4] = j[4].rjust(1) #Astring (see ref link above fn def:)
    for lval in range(len(recarr)):
        lout = [];
        #ATOM/RES details
        lout.append(recarr[lval].ljust(5))#atom#5s [0]
        lout.append(str(aidarr[lval]).rjust(6))#anum#5d/6d [1]
        if len(anamearr[lval]) == 3:
            lout.append(anamearr[lval].rjust(4)) #aname#3s [2]
        else:
            lout.append(anamearr[lval].center(4))#other_anames [2]
        lout.append(resnamearr[lval].ljust(3))#resname#3/4s [3]
        lout.append(str(residarr[lval]).rjust(4))#resid#4s [4] 

        #rxyz details
        if lval != 0 and rxarr[lval] == 0.0 and ryarr[lval] == 0.0\
           and rzarr[lval] == 0.0:
            rx,ry,rz = change_zeros(rxarr[lval-1],ryarr[lval-1],\
                                    rzarr[lval-1],rdist)
            rxarr[lval] = rx; ryarr[lval] = ry; rzarr[lval] = rz

        lout.append(str('%8.3f'%(float(rxarr[lval]))).rjust(8))#x[5]
        lout.append(str('%8.3f'%(float(ryarr[lval]))).rjust(8))#y[6]
        lout.append(str('%8.3f'%(float(rzarr[lval]))).rjust(8))#z[7]
        
        #extra details: occ/beta/elem
        lout.append(str('%6.2f'%(float(occarr[lval]))).ljust(6))#occ[8]
        lout.append(str('%6.2f'%(float(betaarr[lval]))).ljust(6))#temp[9]
        lout.append(segnamearr[lval].ljust(4))#segname[10]
        lout.append(elemarr[lval].rjust(2))#elem[11]

        #write to file
        fout.write("%s%s %s %s %s    %s%s%s%s%s      %s%s\n" \
                   %(lout[0],lout[1],lout[2],lout[3],lout[4],lout[5]\
                     ,lout[6],lout[7],lout[8],lout[9],lout[10],lout[11]))
    fout.close()
#------------------------------------------------------------------
# change zeros to r + rdist*ran
def change_zeros(rxold,ryold,rzold,rmax):
    rval  = random.random()
    theta = 2*math.pi*rval
    phi   = math.pi*rval
    rxnew = rxold + random.random()*rmax*math.sin(phi)*math.cos(theta)
    rynew = ryold + random.random()*rmax*math.sin(phi)*math.sin(theta)
    rznew = rzold + random.random()*rmax*math.cos(phi)
    return rxnew,rynew,rznew
#------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------
