# To generate itp files for certain polymers and acetylated cellulose
# Using Martini V3.0 for cellulose structure
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
    return at_list, bo_list, an_list, di_list
    # Close file
    fpetg.close()
#----------------------------------------------------------------------------------
# Design PLA
def design_pla(nmons,moltype):

    # Ref: MartiniV2.0-Lipids and Li et al., Acta mechanica solida sinica 30 (2017) 630Ã¢€“637
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
        else: #N4a-N4a-N4a (Li et al., Acta mechanica solida sinica 30 (2017) 630-637)
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
               'sinica 30 (2017) 630Ã¢€“637\n')
    
    # Write topology
    write_topology(fpla,moltype,nmons,natoms,nbonds,0,nangls,0,0,\
                   at_list,bo_list,[],an_list,[],[])

    # Close file
    return at_list, bo_list, an_list, di_list
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
    return at_list, bo_list, an_list, di_list
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
#    fin.write("; Include Position restraint file\n")
#    fin.write("#ifdef POSRES\n")
#    fin.write("#include \"posre-%s.itp\"\n" %(moltype.lower()))
#    fin.write("#endif\n")

#----------------------------------------------------------------------------------
