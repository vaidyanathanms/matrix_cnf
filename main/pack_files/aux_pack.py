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
#-------------------------------------------------------------------
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

#Check pdb/psf/top files for the melt (or polymer)
def check_inp_files(dum_inpdir, inp_name):
    # check structure files 
    if glob.glob(dum_inpdir+'/' + inp_name) == []:
        raise RuntimeError(inp_name + " not found in " + dum_inpdir)
#-------------------------------------------------------------------

# Create all output directories
def create_output_dirs(superdir,acetval,acetper,addpoly):

    if not os.path.isdir(superdir):
        os.mkdir(superdir)

    if acetper != 0 and addpoly == 'None':
        outsub = 'acet_cnf_m'+str(acetval)+'_'+str(acetper)
    elif acetper != 0 and addpoly != 'None':
        outsub =  addpoly+'_cnf_m'+ str(acetval)+'_'+str(acetper)
    elif acetper == 0 and addpoly != 'None':
        outsub = addpoly

    subdir = superdir + '/' + outsub
    if not os.path.isdir(subdir):
        os.mkdir(subdir)

    return subdir
#-------------------------------------------------------------------
# Check for gaussian chains and write
def check_gaussianity_and_write(gmxdir,inpfyle2,nmons,nchains,\
                                matname,tol,packalldir):
    

    if not os.path.isdir(packalldir):
        raise RuntimeError(packalldir + ' not found')

    frg = open(packalldir + '/inprgstats_' + matname + '.dat','w')
    frg.write('%s\t%s\t%s\t%s\t%s\t%s\n' 
              %('ChainID','Nres (Nmon)', 'Natoms', \
                'Rg^2','Rg^4','Rg^4/Rg^2'))
    outdir = packalldir + '/' +  'gaussianchains'
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
    write_gaussian_chains(rg2ch,rg4ch,atarr,outdir,inpfyle,matname,tol)
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
def write_gaussian_chains(rg2ch,rg4ch,atarr,dirname,fname,matname,\
                          tolval):
    with open(fname) as fin:
        for ctr in range(len(rg2ch)):
            if abs(1.0-(rg4ch[ctr]/(rg2ch[ctr]**2))) < tolval:
                line = fin.readline()
                with open(dirname+'/'+matname+'_'+str(ctr+1)+'.pdb','w') as fout:
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
    fout.write(' inside box %g  %g  %g  %g  %g  %g\n'
               %(xmin-0.25*dr,ymin-0.25*dr,zmin-dr,xmax+0.25*dr,ymax+0.25*dr,zmax+dr))
    fout.write('end structure\n')
    fout.write('\n')
#------------------------------------------------------------------
# Add polymer matrix
def pack_polymer_matrix(matdir,matname,nch,xmin,ymin,zmin,xmax,ymax,zmax,fout,mag,dmax):
    if not os.path.isdir(matdir):
        raise RuntimeError(matdir + " not found")

    matlist = glob.glob(matdir + '/' + matname + '*.pdb')
    if matlist == []:
        raise RuntimeError('No Gaussian inputs found. Consider changing gaus_tol')
    
    lenlist = len(matlist)
    nsets   = int(nch/lenlist) if int(nch/lenlist) > 1 else 1
    nextra  = nch%lenlist
    dr = (mag-1)*dmax

    for chain in range(min(lenlist,nch)):
        fout.write('structure  %s\n' %(matlist[chain]))
        fout.write('  number   %d\n' %(nsets))  
#        fout.write('  center \n') #Don't center matrix chains
        fout.write(' inside box %g  %g  %g  %g  %g  %g\n'
                   %(xmin-0.5*dr,ymin-0.5*dr,zmin-0.5*dr,xmax+0.5*dr,ymax+0.5*dr,zmax+0.5*dr))
        fout.write('end structure\n')
        fout.write('\n')

    if nextra != 0:
        fout.write('structure  %s\n' %(matlist[0]))
        fout.write('  number   %d\n' %(nextra))  
#        fout.write('  center \n') #Don't center matrix chains
        fout.write(' inside box %g  %g  %g  %g  %g  %g\n'
                   %(xmin-dr,ymin-dr,zmin-dr,xmax+dr,ymax+dr,zmax+dr))
        fout.write('end structure\n')
        fout.write('\n')
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
def make_fin_dirs_and_copy(scr_dir,matrix,mod_cell,run_pack,packfyle):
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
        raise RuntimeError('No patch for acetylation length: ', acet_val)
#------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------
