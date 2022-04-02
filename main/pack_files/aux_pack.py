# Auxiliary file for create_packmol_file
import os
import sys
import numpy as np
import re
import shutil
import glob
import math
import struct
import subprocess
import fileinput
#-------------------------------------------------------------------
# Generic copying script
def gencpy(dum_maindir,dum_destdir,fylname):
    
    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        raise RuntimeError('ERROR: ', srcfyl, 'not found')

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#-------------------------------------------------------------------
# Check directories
def check_dir(dirname):
    if not os.path.isdir(dirname): # 
        print("FATAL ERROR: ", dirname, " not found\n")
        raise RuntimeError("Check directory path: " + dirname)
#-------------------------------------------------------------------
#Find pdb file for the matrix/additive
def find_inp_files(dum_inpdir, poltype, maindir):
    os.chdir(dum_inpdir)
    if glob.glob('*.pdb') == []:
        raise RuntimeError("No pdb files formed found for "+poltype)
    elif len(glob.glob('*.pdb')) == 1:
        pdb_fname = glob.glob('*.pdb')
    else:
        print("Multiple PDB files found. Using latest for "+poltype)
        fnames = glob.glob('*.pdb')
        pdb_fname = max(fnames, key=os.path.getctime)
    os.chdir(maindir)
    return pdb_fname[0]
#-------------------------------------------------------------------
# Create all output directories
def create_output_dirs(superdir,acetval,acetper,addpoly):
    if not os.path.isdir(superdir):
        os.mkdir(superdir)

    if acetper != 0 and addpoly.lower() == 'None'.lower():
        outsub = 'acet_cnf_m'+str(acetval)+'_'+str(acetper)
    elif acetper != 0 and addpoly.lower() != 'None'.lower():
        outsub = 'acet_cnf_m'+ str(acetval)+'_'+str(acetper)\
                 + '_with_' + addpoly
    elif acetper == 0 and addpoly.lower() != 'None'.lower():
        outsub = 'with_' + addpoly

    subdir = superdir + '/' + outsub
    if os.path.isdir(subdir):
        clear_all_files(subdir)    
    if not os.path.isdir(subdir):  
        os.mkdir(subdir)
    return subdir
#-------------------------------------------------------------------
# Clear all files from the directory
def clear_all_files(subdir):
    print("Output directory already exists.")
    inp = input("Clean and rewrite directory (y/n)? ")
    while not inp.lower() == 'y' or inp.lower() == 'n':
        print("Please input y or n: ")
        inp = input("Clean directory and rewrite files (y/n)? ")
    if inp.lower() == 'y':
        shutil.rmtree(subdir)
    else:
        quit()
#-------------------------------------------------------------------
# Check for gaussian chains and write
def check_gaussianity_and_write(gmxdir,inpfyle2,nmons,nchains,\
                                polyname,tol,packalldir):
    

    if not os.path.isdir(packalldir):
        raise RuntimeError(packalldir + ' not found')

    frg = open(packalldir + '/inprgstats_' + polyname + '.dat','w')
    frg.write('%s\t%s\t%s\t%s\t%s\t%s\n' 
              %('ChainID','Nres (Nmon)', 'Natoms', \
                'Rg^2','Rg^4','Rg^4/Rg^2'))
    outdir = packalldir + '/' + polyname + '_gaussianchains'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        for f in glob.glob(outdir + '/*'):
            os.remove(f)

    inpfyle = gmxdir + '/' + inpfyle2

    with open(inpfyle) as fin:
        rx = []; ry = []; rz = []
        rg2ch = []; rg4ch = []; atarr = []
        chcntr = 1; atcnt = 0
        for line in fin:
            if line.lstrip().startswith('#') or \
               line.lstrip().startswith('@'):
                continue
            elif 'TER' in line:
                if mval != nmons:
                    raise RuntimeError("# of monomers (residues) mismatch",\
                                       chcntr,mval,nmons)
                rg2, rg4 = compute_rg(rx,ry,rz,atcnt)
                rg2ch.append(rg2); rg4ch.append(rg4)
                atarr.append(atcnt)
                frg.write('%g\t%g\t%g\t%g\t%g\t%g\n' \
                          %(chcntr,mval,atcnt,rg2,rg4,rg4/rg2**2))
                rx = []; ry = []; rz = []
                chcntr += 1; atcnt = 0
            else:
                mval = int(line.split()[4])
                rx.append(float(line.split()[5]))
                ry.append(float(line.split()[6]))
                rz.append(float(line.split()[7]))
                atcnt += 1
    write_gaussian_chains(rg2ch,rg4ch,atarr,outdir,inpfyle,polyname,tol)
    print('Finished guassianity writes...')
    rg2max = max(rg2ch)
    return outdir,math.sqrt(rg2max)
#------------------------------------------------------------------
# Compute Rg
def compute_rg(rx, ry, rz, atcnt):

    rxa = np.array(rx); rya = np.array(ry); rza = np.array(rz)
    rxcm = np.average(rxa); rycm = np.average(rya)
    rzcm = np.average(rza)

    rx2 = np.sum(np.square(rxa-rxcm)); rx4 = np.sum(np.power(rxa-rxcm,4)); 
    ry2 = np.sum(np.square(rya-rycm)); ry4 = np.sum(np.power(rya-rycm,4)); 
    rz2 = np.sum(np.square(rza-rzcm)); rz4 = np.sum(np.power(rza-rzcm,4)); 

    rx2 /= atcnt; ry2 /= atcnt; rz2 /= atcnt
    rx4 /= atcnt; ry4 /= atcnt; rz4 /= atcnt

    return rx2 + ry2 + rz2, rx4 + ry4 + rz4
#------------------------------------------------------------------
# Write Gaussian Chains
def write_gaussian_chains(rg2ch,rg4ch,atarr,dirname,fname,polyname,\
                          tolval):
    with open(fname) as fin:
        for ctr in range(len(rg2ch)):
            if abs(1.0-(rg4ch[ctr]/(rg2ch[ctr]**2))) < tolval:
                line = fin.readline()
                with open(dirname+'/'+polyname+'_'+str(ctr+1)+'.pdb','w') as fout:
                    fout.write('REMARK 999 chainID: %d; Rg4/(Rg2)^2: %g\n' \
                               %(ctr+1,rg4ch[ctr]/(rg2ch[ctr]**2)))
                    fout.write(line)
                    while not 'TER' in line:
                        line = fin.readline()
                        fout.write(line)
            else:
                while not 'TER' in fin.readline():
                    pass
#------------------------------------------------------------------
# Make packmol headers
def packmol_headers(fout,matname,outdir,acetper,acetval,addpoly):
    outname = 'packed_'+matname+'.pdb'
    fout.write('# Pack cellulose and polymer matrix\n')
    fout.write('# Inputs\n\n')
    fout.write('tolerance  %g\n' %(2.0))
    if acetper != 0 and addpoly == 'None':
        outfyle = 'packed_'+matname+'_cnf_m'+str(acetval)+'.pdb'
    elif acetper != 0 and addpoly != 'None':
        outfyle = 'packed_'+matname+'_'+addpoly+'_cnf_m'+\
                  str(acetval)+'.pdb'
    elif acetper == 0 and addpoly != 'None':
        outfyle = 'packed_'+matname+'_'+addpoly+'.pdb'
        
    fout.write('filetype  pdb\n')
    fout.write('output  %s\n' %(outdir + '/' + outfyle))
    fout.write('seed  %d\n' %(-1))
    fout.write('\n#Begin structure generation \n\n')
    return outname
#------------------------------------------------------------------
# Find cellulose/modified cellulose dimensions
def read_pdb_file(inpfyle):
    fmt_at = b'4sx6sx4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s'
    fmt_atq= b'4sx6sx4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s2s'
    fmt_ht = b'6s5sx4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s'
    fmt_htq= b'6s5sx4s1s4s1s4s1s3x8s8s8s6s6s6x4s2s2s'
    aidarr = []; anamearr = []; residarr = []
    rxarr = []; ryarr = []; rzarr = []
    segnamearr = [];massarr = []
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
                    = struct.unpack(fmt_atq,line.strip())
            elif line[0:4] == 'HEATOM' and len(line.strip()) == 78:
                (rec,aid,aname,altloc,resname,chid,resid,rescode,\
                 rx,ry,rz,occ,beta,segname,elem) \
                    = struct.unpack(fmt_ht,line.strip())
            elif line[0:4] == 'HEATOM' and len(line.strip()) == 80:
                (rec,aid,aname,altloc,resname,chid,resid,rescode,\
                 rx,ry,rz,occ,beta,segname,elem,qval) \
                    = struct.unpack(fmt_htq,line.strip())
            elif line[0:3] == 'END' or line[0:3] == 'TER':
                continue
            else:
                raise RuntimeError("Unknown line format", line)
            
            aidarr.append(int(aid.decode('utf-8')))
            anamearr.append(aname.decode('utf-8'))
            residarr.append(int(resid.decode('utf-8')))
            rxarr.append(float(rx.decode('utf-8')))
            ryarr.append(float(ry.decode('utf-8')))
            rzarr.append(float(rz.decode('utf-8')))
            segnamearr.append(segname.decode('utf-8'))
            
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
            massarr.append(mval)

    return aidarr,residarr,rxarr,ryarr,rzarr,segnamearr,massarr
#------------------------------------------------------------------
# Return cnf max and min dimensions
def find_cnf_dim(acetdir,acetfyle2,packdir):   
    gencpy(packdir,acetdir,acetfyle2+'.pdb')
    gencpy(packdir,acetdir,acetfyle2+'.psf')
    acetfyle = acetdir + '/' + acetfyle2 + '.pdb'
    aidcnf,residcnf,rxcnf,rycnf,rzcnf,segnamecnf,masscnf \
        = read_pdb_file(acetfyle)
#    rxcm,rycm,rzcm = find_com_chains(residcnf,masscnf,rxcnf,rycnf,rzcnf)
    return min(rxcnf),min(rycnf),min(rzcnf),max(rxcnf),max(rycnf),\
        max(rzcnf)
#------------------------------------------------------------------
# return com of chains
def find_r_com_chains(resarr,massarr,rxarr,ryarr,rzarr):
    chid = 1;ctr = 1
    rxcmarr = []; rycmarr = []; rzcmarr = []
    rxcm = 0.0; rycm = 0.0; rzcm = 0.0
    while chid <= max(resarr):
        if resarr[ctr] == chid:
            rxcm += massarr[ctr]*rxarr[ctr]
            rycm += massarr[ctr]*ryarr[ctr]
            rzcm += massarr[ctr]*rzarr[ctr]
            ctr += 1
        else:
            chid += 1
            rxcmarr.append(rxcm)
            rycmarr.append(rycm)
            rzcmarr.append(rzcm)
            rxcm = 0.0; rycm = 0.0; rzcm = 0.0
    return rxcmarr, rycmarr, rzcmarr
#------------------------------------------------------------------
# Add cellulose chains
def pack_cellulose_chains(fout,packdir,inpfyle,ncnf,xmin,ymin,zmin,xmax,ymax,zmax,mag,dmax):
    dr = (mag-1)*dmax
    fout.write('structure  %s\n' %(packdir+'/'+inpfyle+'.pdb'))
    fout.write('  number   %d\n' %(ncnf))  
    fout.write('  center \n')
#    fout.write(' inside box %g  %g  %g  %g  %g  %g\n'
#               %(xmin-dr,ymin-dr,zmin-dr,xmax+dr,ymax+dr,zmax+dr))
    fout.write(' fixed %g  %g  %g  %g  %g  %g\n'
               %(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)) #center to zero
    fout.write('end structure\n')
    fout.write('\n')
#------------------------------------------------------------------
# Add main polymer matrix/additional polymer blend/co/homoplymers
def pack_polymers(matdir,matname,nch,xmin,ymin,zmin,xmax,ymax,\
                  zmax,fout,mag,dmax):
    if not os.path.isdir(matdir):
        raise RuntimeError(matdir + " not found")

    matlist = glob.glob(matdir + '/' + matname + '*.pdb')

    if matlist == []:
        raise RuntimeError('No Gaussian inputs found. Consider changing gaus_tol')
    
    lenlist = len(matlist)
    if nch/lenlist > 1:
        nsets = int(nch/lenlist); nextra = nch%lenlist
    else:
        nsets = 1; nextra = 0
    dr = (mag-1)*dmax
    
    lxh = (xmax-xmin)/2; lyh = (ymax-ymin)/2; lzh = (zmax-zmin)/2
    for chain in range(min(lenlist,nch)):
        fout.write('structure  %s\n' %(matlist[chain]))
        fout.write('  number   %d\n' %(nsets))  
#        fout.write('  center \n') #Don't center matrix chains
        fout.write(' inside box %g  %g  %g  %g  %g  %g\n'\
                   %(-lxh-dr,-lyh-dr,-lzh-dr,lxh+dr,lyh+dr,lzh+dr))
#                   %(xmin-dr,ymin-dr,zmin-dr,xmax+dr,ymax+dr,zmax+dr))
        fout.write('end structure\n')
        fout.write('\n')

    if nextra != 0:
        fout.write('structure  %s\n' %(matlist[0]))
        fout.write('  number   %d\n' %(nextra))  
#        fout.write('  center \n') #Don't center matrix chains
        fout.write(' inside box %g  %g  %g  %g  %g  %g\n'\
                   %(-lxh-dr,-lyh-dr,-lzh-dr,lxh+dr,lyh+dr,lzh+dr))
        fout.write('end structure\n')
        fout.write('\n')
#------------------------------------------------------------------
# Consistency checks for additional polymers
def expoly_check(ex_ptype,ex_nch,ex_nmon,addpoly,ex_gaus):
    if addpoly.lower() == 'blend':
        if len(ex_ptype) != 2:
            raise RuntimeError('Blends should have exactly two poly types')
    elif addpoly.lower() == 'homo':
        if len(ex_ptype) != 1:
            raise RuntimeError('Homopolymers can have only one poly type')
    elif addpoly.lower() == 'copoly':
        if len(ex_ptype) != 0:
            raise RuntimeError('Copoly should have at least one poly type')
    else:
        raise RuntimeError('Unknown add_poly option: '+ addpoly)

    if len(ex_ptype) != len(ex_nch) or len(ex_ptype) != len(ex_nmon)\
       or len(ex_ptype) != len(ex_gaus):
        raise RuntimeError('Size mismatch between polymer types ' +
                           'and number of chains/monomers/tol')
#------------------------------------------------------------------
# Run packmol
def run_packmol(packfyle,pack_exe,destdir,packsh,currdir):
    if not os.path.exists(pack_exe):
        raise RuntimeError("Packmol executable not found\n")

    print("Copying packmol run script..")
    if not os.path.exists(currdir + '/' + packsh):
        raise RuntimeError(packsh + " not found in " + currdir)
    else:
        gencpy(currdir,destdir,packsh)
       
    os.chdir(destdir)
    print("Running packmol using", packfyle)
    
    rev_fname = packsh.replace('_pyinp','')
    fr  = open(destdir + '/' + packsh,'r')
    fw  = open(destdir + '/' + rev_fname,'w')
    fid = fr.read().replace("py_packexe",pack_exe).\
          replace("py_packinp",packfyle)

    fw.write(fid)
    fw.close()
    fr.close()
    
    if not os.path.isdir("outdir"):
        os.mkdir("outdir")

    subprocess.call(["sbatch", rev_fname])
#------------------------------------------------------------------
# Make final directories
def make_fin_dirs_and_copy(scr_dir,matrix,mod_cell,run_pack,\
                           packfyle):
    if not os.path.isdir(scr_dir):
        os.mkdir(scr_dir)
    poly_dir = scr_dir + '/' + matrix
    if not os.path.isdir(poly_dir):
        os.mkdir(poly_dir)
    if mod_cell:
        poly_dir = poly_dir + '/mod_cell'
        if not os.path.isdir(poly_dir):
            os.mkdir(poly_dir)
    if run_pack:
        if not os.path.exists(packfyle):
            raise RuntimeError(packfyle + " not found\n")
        gencpy(os.getcwd(),poly_dir,packfyle)
#------------------------------------------------------------------
# Make acetylated cellulose
def make_acet_cell(acet_dir,acet_val,cell_dp,acet_per,acet_tol,\
                   acet_fyle,acet_att,pack_mat,acet_new,nativ_cnf):
    if acet_new:
        if not os.path.exists(acet_dir + '/genmodify_cnf_pyinp.tcl'):
            print('ERR: genmodify_cnf_pyinp.tcl not found\n')
            raise RuntimeError('File to generate modified cellulose not found\n')


        if glob.glob(acet_fyle+'*'):
            for fyle in glob.glob(acet_fyle):
                os.remove(fyle)

        acet_pat  = ret_acet_pat(acet_val)
        tcl_fname = 'genmodify_cnf_pyinp.tcl'
        rev_fname = tcl_fname.replace('_pyinp','')
        fr  = open(acet_dir + '/' + tcl_fname,'r')
        fw  = open(pack_mat + '/' + rev_fname,'w')
        fid = fr.read().replace("py_nmons",str(cell_dp)).\
              replace("py_acprob",str(acet_per)).\
              replace("py_relerr",str(acet_tol)).\
              replace("py_maxatt",str(acet_att)).\
              replace("py_outname",str(pack_mat + '/' + acet_fyle)).\
              replace("py_patch",str(acet_pat)).\
              replace("py_acetval",str(acet_val)).\
              replace("py_outdir",str(pack_mat)).\
              replace("py_cnfmain",str(nativ_cnf))
        fw.write(fid)
        fw.close()
        fr.close()

        os.chdir(pack_mat)
        subprocess.call(['vmd','-dispdev','text','-e',rev_fname])

    # copy acetylated cellulose to final directory
    lst_fyle = glob.glob(pack_mat + '/' + acet_fyle+'*')
    if len(lst_fyle) == 0:
        raise RuntimeError('No acetylated chains produced\n')
#------------------------------------------------------------------
# Return name of acetylation patch
def ret_acet_pat(acet_val):
    if acet_val == 1:
        return '6TAC'
    elif acet_val == 3:
        return '6AC3'
    elif acet_val == 7:
        return '6AC7'
    elif acet_val == 11:
        return '6AC11'
    else:
        raise RuntimeError('No patch for acetylation length: '+acet_val)
#------------------------------------------------------------------
# Make gmx top file for acetylated cellulose
def make_top_file_for_acetcell(maindir,destdir,topdir,topout):
    # Replace tcl files with relevant inputs
    if not os.path.exists(maindir + '/gen_celltop_pyinp.tcl'):
        raise RuntimeError('gen_celltop_pyinp.tcl not found in '\
                           +maindir)

    fr  = open(maindir + '/gen_celltop_pyinp.tcl','r')
    fw  = open(destdir + '/gen_celltop.tcl','w')
    fid = fr.read().replace("py_topdir",topdir).\
          replace("py_outname",topout)
    
    fw.write(fid)
    fw.close()
    fr.close()
    # Generate topology file
    os.chdir(destdir)
    subprocess.call(['vmd','-dispdev','text','-e','gen_celltop.tcl'])
    if not os.path.exists(destdir + '/' + topout + '.top'):
        raise RuntimeError(topout + '.top not found in ' + destdir)
    os.chdir(maindir)
#------------------------------------------------------------------
# Split gmx top file into prm and itp files
def split_top_file_to_prmitp(top_file,destdir,maindir):
    os.chdir(destdir)
    mol_info = []
    if not os.path.isdir(destdir + '/all_toppar'):
        os.mkdir(destdir + '/all_toppar')
    prm_file = destdir + '/all_toppar/' + top_file + '.prm'
    itp_file = destdir + '/all_toppar/' + top_file + '.itp'
    # Write until [ moleculetype ] to *.prm file and
    # Write from [ moleculetype ] until [ System ] to *.itp
    # Then skip from [ System ] to [ molecules ]
    # Finally add system info to arrays
    with open(top_file+'.top','r') as ftop, \
         open(prm_file,'w') as fprm,\
         open(itp_file,'w') as fitp:
        line = ftop.readline()
        while not '[ moleculetype ]' in line.lower():
            fprm.write(line)
            line = ftop.readline()
        while not '[ system ]' in line.lower():
            fitp.write(line)
            line = ftop.readline()
        while not '[ molecules ]' in line.lower():
            line = ftop.readline()
        mol_info.append(';cellulose_main \n')
        for line in ftop:
            if not line.lstrip().startswith(';'):
                mol_info.append(line)
    os.chdir(maindir)
    return prm_file, itp_file, mol_info
#------------------------------------------------------------------
# Write headers to alltop.top
def write_header_topfile(f_all):
    f_all.write(';; -------------------------------------\n')
    f_all.write(';; Combined topology file\n')
    f_all.write(';; Generated using create_packmol_file.py\n')
    f_all.write(';; Use with GROMACS grompp \n \n')
    f_all.write(';; -------------------------------------\n\n')
#------------------------------------------------------------------
# Copy forcefield/top files of matrix to packmol directory
#polygmxdir - original CHARMM-GUI directory
#pack_dir   - combined output directory
def copy_pol_toppar_files(polgmxdir,pack_dir,prmfyle_arr,\
                          itpfyle_arr,mol_info,mol_infoptr,moltyp):

    maindir = os.getcwd()
    # First copy files into pack_dir
    if not os.path.isdir(polgmxdir+'/toppar'):
        raise RuntimeError('toppar directory not found in '\
                           + polgmxdir)
    if not os.path.isdir(pack_dir + '/all_toppar'):
        os.mkdir(pack_dir + '/all_toppar')

    os.chdir(polgmxdir + '/toppar')
    mat_toppar = glob.glob('*')
    if mat_toppar == []:
        raise RuntimeError('No topology files found in '\
                           + polgmxdir + '/toppar')

    # Change file name so that files are not overwritten
    for fyle in mat_toppar:
        srcfyle = polgmxdir + '/toppar/' + fyle
        desfyle = pack_dir + '/all_toppar/' + moltyp + fyle
        shutil.copy2(srcfyle,desfyle)

    # Now all files are in pack_dir + '/all_toppar'
    # Find prm and itp files from this folder
    # File names are changed with the moltype prefix
    # This way original files are intact
    # Read toppar files from pack_dir + '/all_toppar'

    os.chdir(pack_dir + '/all_toppar')
    itp = -1; prm = -1
    for fyle in mat_toppar:
        fname = moltyp + fyle
        if not os.path.exists(fname):
            raise RuntimeError(fname+' not found in ',os.getcwd())
        with open(fname,'r') as fin: # check filetype is itp or prm
            for line in fin:
                if '[ moleculetype ]' in line:
                    itpfyle_arr.append(pack_dir +\
                                       '/all_toppar/' + fname)
                    itp = 1
                    break
                elif '[ atomtypes ]' in line or \
                     '[ bondtypes ]' in line:
                    prmfyle_arr.append(pack_dir + \
                                       '/all_toppar/' + fname)
                    prm = 1
                    break

    # Sanity check to see whether itp/prm files are found
    # Checking fname or fyle above is equivalent because
    # if fyle is found in polygmxdir/toppar, 
    # then fname should be present in pack_dir/all_toppar

    if itp == -1:
        raise RuntimeError('Topology files not found in '\
                           + polgmxdir + '/toppar')
    if prm == -1:
        raise RuntimeError('Interaction files not found in '\
                           + polgmxdir + '/toppar')

    # Open topol file to obtain molecule name and size
    # IMPORTANT: However change the moleculename so that multiple
    # molecules with same names are not there.
    # No need to copy this to pack_dir
    if not os.path.exists(polgmxdir + '/topol.top'):
        raise RuntimeError('topol.top missing in ' + polgmxdir)
    else:
        with open(polgmxdir + '/topol.top','r') as ftop:
            while not '[ molecules ]' in ftop.readline():
                pass
            line = ftop.readline()
            mol_info.append(';'+moltyp)#to know new molecule is added
            mol_infoptr.append(len(mol_info))#points to index of molname
            for line in ftop:
                if not line.lstrip().startswith(';'):
                    # set default number of chains as 0. 
                    # set correct value in combine_top_files
                    molname = moltyp + '_' + line.lstrip().split()[0]
                    mol_info.append(molname + '\t' + str(0))
        gencpy(polgmxdir,pack_dir,'topol.top')
    os.chdir(maindir)
#------------------------------------------------------------------
# Change forcefield files to add ;
def add_comment_to_ff_files(prmfyles):
    # Do not edit the first file which is CNF
    for i in range(1,len(prmfyles)): # don't edit first file
        with fileinput.FileInput(prmfyles[i],inplace=True,\
                                 backup='.bak') as fin:
            for line in fin:
                if 'defaults' in line:
                    print('; '+ line,end='')
                    for j in range(2):
                        line = fin.readline()
                        print('; '+ line,end='')
                else:
                    print(line,end='')
#------------------------------------------------------------------
# Change molecule name in the itp (topology) file
def change_molname_itp_file(molinfo,molptr,itpfyles):
    # No need to edit the first file which is CNF 
    for i in range(1,len(itpfyles)): # don't edit first file
        with fileinput.FileInput(itpfyles[i],inplace=True,\
                                 backup='.bak') as fin:
            for line in fin:
                if '[ moleculetype ]' in line:
                    print(line,end='')
                    line = fin.readline()
                    # Read until the next uncommented line
                    while line.lstrip().startswith(';'):
                        print(line,end='')
                        line = fin.readline()
                    # Replace molecule name
                    # Check whether next main line has 2 entries
                    lsplit = line.split() #don't strip \n
                    if len(lsplit) != 2:
                        raise RuntimeError('[ moleculetype ] should contain'\
                                           +' only 2 args in the itp file: '\
                                           + itpfyles[i], lsplit)
                    moltype = molinfo[molptr[i-1]].lstrip().split()[0]
                    lsplit[0] = moltype
                    print('\t'.join(lsplit))
                else:
                    print(line,end='')
#------------------------------------------------------------------
# Combine polymer matrix and acetylation files
def combine_top_files(outdir,prm_files,itp_files,molinfo,nch,\
                      addpoly,exnch):
    all_top  = outdir + '/alltop.top'
    f_all = open(all_top,'w')
    write_header_topfile(f_all)

    # Write Interaction parameter files
    f_all.write(';; Include interaction parameter files\n')
    for i in range(len(prm_files)):
        f_all.write('#include "%s"\n' %(prm_files[i]))
    f_all.write('\n')

    # Write Topology parameter files
    f_all.write(';; Include topology data files\n')
    for i in range(len(itp_files)):
        f_all.write('#include "%s"\n' %(itp_files[i]))
    
    # Write System data
    f_all.write('\n')
    f_all.write('[ System ]\n')
    f_all.write('; Name\n')
    f_all.write(' Combined\n')

    # Write molecules data
    f_all.write('\n')
    f_all.write('[ molecules ]\n')
    f_all.write('; %s\t%s\n' %('Compound', '#mols'))

    # Write cellulose part
    i = 0
    while ';matrix' not in molinfo[i] and i < len(molinfo):
        f_all.write('%s' %(molinfo[i])); i+=1
            
    # Write matrix part
    while i < len(molinfo):
        if addpoly in molinfo[i]:
            break
        if ';matrix' in molinfo[i]:
            f_all.write('%s\n' %(molinfo[i])); i+=1
        else:
            moltype = molinfo[i].lstrip().split()[0]
            f_all.write('%s\t%d' %(moltype,nch)); i+=1

    # Write additional polymers
    if addpoly.lower() != 'None'.lower():
        j = 0
        f_all.write('\n') # one line is necessary
        while i < len(molinfo):
            if ';'+addpoly in molinfo[i]:
                f_all.write('%s\n' %(molinfo[i])); i+=1
            else:
                moltype = molinfo[i].lstrip().split()[0]
                f_all.write('%s\t%d\n' %(moltype,exnch[j]));
                i+=1; j+=1

    f_all.close()
    return all_top
#------------------------------------------------------------------
# Clean and sort files
def clean_and_sort_files(outdir,acetpref):
    maindir = os.getcwd()
    move_ftype('pdb',outdir,acetpref,maindir)
    move_ftype('psf',outdir,acetpref,maindir)
    move_ftype('tcl',outdir,acetpref,maindir)
    move_ftype('dat',outdir,acetpref,maindir)
#------------------------------------------------------------------
# To move file of a particular type
def move_ftype(ftype,outdir,acetpref,currdir):
    os.chdir(outdir)
    dname = outdir + '/' + ftype + '_files'
    if not os.path.isdir(dname):
        os.mkdir(dname)
    for fyle in glob.glob('*' + ftype):
        if 'packed' in fyle or acetpref in fyle:
            continue
        else:
            shutil.move(fyle,dname+'/' + fyle)
#------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------
