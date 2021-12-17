# Auxiliary file to manipulate and edit GROMACS files
#------------------------------------------------------------------

#Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
#------------------------------------------------------------------
# Set default thermostat coefficients
def couple_coeff(inp_type,coeff_fyle = 'None'):
    # default. change if needed
    ref_pres = 1
    tau_temp_nvt     = 0.1
    tau_temp_nvthigh = 0.2
    tau_temp_berend  = 0.1
    tau_temp_parrah  = 0.2
    tau_pres_berend  = 0.5
    melt_topname     = 'None'
    if inp_type == 'melts':
        tau_pres_parrah  = 5.0
    else:
        tau_pres_parrah  = 8.0
    if coeff_fyle != 'None':
        with open(coeff_fyle) as farg:
            for line in farg:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    continue
                words = line.split()
                if words[0] == 'TauTemp_NVT':
                    tau_temp_nvt  = float(words[1])
                elif words[0] == 'TauTempHigh_NVT':
                    tau_temp_nvthigh  = float(words[1])
                elif words[0] == 'TauTemp_Berendsen':
                    tau_temp_berend = float(words[1])
                elif words[0] == 'TauPres_Berendsen':
                    tau_pres_berend = float(words[1])
                elif words[0] == 'TauTemp_Parrinello':
                    tau_temp_parrah = float(words[1])
                elif words[0] == 'TauPres_Parrinello':
                    tau_pres_parrah = float(words[1])
                elif words[0] == 'Ref_Pres':
                    ref_pres  = float(words[1])
                elif words[0] == 'Melt_Topfile':
                    melt_topname = words[1]
                else:
                    raise RuntimeError("Unknown keyword: "+ words[0] \
                                       + " in " + str(coeff_fyle))
    return tau_temp_nvt,tau_temp_nvthigh,tau_temp_berend,\
        tau_temp_parrah,tau_pres_berend,tau_pres_parrah,ref_pres,\
        melt_topname
#------------------------------------------------------------------

#Check pdb/psf/top files for the melt (or polymer)
def check_inp_files(dum_inpdir,top_name):
    # check structure files (*.pdb/.gro)
    if glob.glob(dum_inpdir+'/*.pdb') == [] and \
       glob.glob(dum_inpdir+'/*.gro') == []:
        raise RuntimeError("No polymer pdb/gro files found in "\
                           + str(dum_inpdir))
    elif len(glob.glob(dum_inpdir+'/*.gro')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.gro')[0]
    elif len(glob.glob(dum_inpdir+'/*.pdb')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.pdb')[0]
    elif len(glob.glob(dum_inpdir+'/*.gro')) > 1:
        print('More than one config file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.gro')
        conf_fname = max(fnames, key=os.path.getctime)
    elif len(glob.glob(dum_inpdir+'/*.pdb')) > 1:
        print('More than one config file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.pdb')
        conf_fname = max(fnames, key=os.path.getctime)

    # check topology files
    if glob.glob(dum_inpdir+'/*.top') == []:
        raise RuntimeError("No polymer topology files found in "\
                           +str(dum_inpdir))
    if top_name != 'None':
        if not os.path.exists(dum_inpdir + '/' + top_name):
            raise RuntimeError("Specified poly top file not found")
        else:
            topol_fname = top_name
    elif len(glob.glob(dum_inpdir+'/*.top')) == 1 and \
         top_name == 'None':
        topol_fname = glob.glob(dum_inpdir+'/*.top')[0]
    else:
        print('More than one topology file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.top')
        topol_fname = max(fnames, key=os.path.getctime)
    
    return conf_fname, topol_fname
#------------------------------------------------------------------

# Find high temperature equilibrated file for sub-high_T systems
def find_hightemp_conf(workdir1,high_temp):

    hi_T_dir = workdir1 + '/T_'+str(high_temp)
    if not os.path.isdir(hi_T_dir):
        print("WARNING: Could not find high temperature directory at T= " + \
              str(high_temp) + "\n This is required for low temperature simulations." + \
              "\n Using default configuration.\n")
        return 'none'

    hiT_conf_files = glob.glob(hi_T_dir + '/*.gro')
    if hiT_conf_files == []:
        print("WARNING: Could not find any gro files for high temperature; T= " + \
              str(high_temp) + "\n This is required for low temperature simulations." + \
              "\n Using default configuration.\n")
        return 'none'
        
    conf_fname = max(hiT_conf_files, key=os.path.getctime)
    poly_conf = ret_file_str(conf_fname)
    gencpy(hi_T_dir,workdir1,poly_conf) #copy to head directory
    return conf_fname #return full path
#------------------------------------------------------------------

# Copy polyconf/polytop for this for this temp from main workdir
def cpy_polyfiles_to_tempdir(poly_conffile,poly_topfile,workdir1,\
                             temp_dir):
    poly_conf = ret_file_str(poly_conffile)
    gencpy(workdir1,temp_dir,poly_conf)
    poly_top = ret_file_str(poly_topfile)
    gencpy(workdir1,temp_dir,poly_top)
#------------------------------------------------------------------

# Create index groups
def create_indxgrps(destdir,inp_type,npoly_res,solv_name,wat_name):
    tcgrp_fname = 'tcgrp_inp.txt'
    if inp_type == 'melts':
        tc_str = 'System'
    elif inp_type == 'solvents':
        tc_str = 'resnr 1 to '+ str(npoly_res)
        tc_str = tc_str + '; ' + 'resname ' + solv_name + '\n'
    elif inp_type == 'cosolvents':
        tc_str = 'resnr 1 to '+ str(npoly_res)
        tc_str = tc_str + '; ' + 'resname ' + solv_name
        tc_str = tc_str + '; ' + 'resname ' + \
                 wat_name.split('_')[0] + '\n'
    
    with open(destdir+'/'+tcgrp_fname,'w') as fw:
        fw.write(tc_str)
    return tcgrp_fname
#------------------------------------------------------------------

# Create temperature coupling groups
def create_tcgrps(inp_type,npoly_res,solv_name,wat_name):
    if inp_type == 'melts':
        tc_grp = 'System'
        tc_typ = 'Single'
    elif inp_type == 'solvents':
        tc_grp = 'resnr_1_to_'+ str(npoly_res)
        tc_grp = tc_grp + '  ' + 'resname_' + solv_name
        tc_typ = 'Multi '
    elif inp_type == 'cosolvents':
        tc_grp = 'resnr_1_to_'+ str(npoly_res)
        tc_grp = tc_grp + '  ' + 'resname_' + solv_name
        tc_grp = tc_grp + '  ' + 'resname_' + wat_name
        tc_typ = 'Multi '
    return tc_grp, tc_typ
#------------------------------------------------------------------

# Check for mdp files and copy/edit if not present
def check_cpy_mdp_files(srcdir,destdir,mdp_fyles,inp_type,Tetau_nvt\
                        ,Tetau_highnvt,Tetau_berend,Tetau_parrah,\
                        Prtau_berend,Prtau_parrah,ref_temp,hi_ref_t,\
                        ref_pres,tc_grpdata,tc_grptype,headdir,coeff_fyle):

    if coeff_fyle != 'None': #gmx inp file from sys.argv
        print("copying", coeff_fyle)
        gencpy(headdir,destdir,coeff_fyle)

    # Only temperatures need to have multiple coupling groups
    str_Tetau_nvt   = genstr(inp_type,Tetau_nvt)
    str_Tetau_hinvt = genstr(inp_type,Tetau_highnvt)
    str_Tetau_ber   = genstr(inp_type,Tetau_berend)
    str_Tetau_par   = genstr(inp_type,Tetau_parrah)
    str_temp        = genstr(inp_type,ref_temp)
    str_hi_temp     = genstr(inp_type,hi_ref_t)

    for mdp_fname in mdp_fyles:
        if not os.path.exists(srcdir + '/' + mdp_fname):
            raise RuntimeError(mdp_fname," not found in ",srcdir)
        py_fname = mdp_fname
        rev_fname = mdp_fname.replace('_pyinp','')
        fr  = open(srcdir + '/' + py_fname,'r')
        fw  = open(destdir + '/' + rev_fname,'w')
        fid = fr.read().replace("py_tcgrps",tc_grpdata).\
              replace("py_grptype",tc_grptype).\
              replace("py_Temptau_vr",str_Tetau_nvt).\
              replace("py_HighTemptau_vr",str_Tetau_hinvt).\
              replace("py_Temptau_Berend", str_Tetau_ber).\
              replace("py_Temptau_ParRah",str_Tetau_par).\
              replace("py_Prestau_Berend",str(Prtau_berend)).\
              replace("py_Prestau_ParRah",str(Prtau_parrah)).\
              replace("py_ref_t",str_temp).\
              replace("py_Highref_t",str_hi_temp).\
              replace("py_ref_p",str(ref_pres))
        fw.write(fid)
        fw.close()
        fr.close()
#------------------------------------------------------------------

# Generate strings
def genstr(inp_type,inp_vals):
    if inp_type == 'melts':
        out_str = str(inp_vals)
    elif inp_type == 'solvents':
        out_str = str(inp_vals) + '  ' + str(inp_vals)
    elif inp_type == 'cosolvents':
        out_str = str(inp_vals) + '  ' + str(inp_vals) + '  ' + \
                  str(inp_vals)
    return out_str
#------------------------------------------------------------------

# Copy shell script files
def cpy_sh_files(srcdir,destdir,sh_pp_fyle,sh_md_fyle):

    # Check for tpr files to provide run conditions
    if glob.glob(destdir+'/*.tpr') == []:
        print('No tpr files found. Beginning new runs..')
        continue_run = 0
    else:
        print('tpr files found. Continuing runs..')
        continue_run = 1

    if continue_run == 0: 
        if not os.path.exists(srcdir+'/' + sh_pp_fyle):
            raise RuntimeError(sh_pp_fyle+" not found in "+\
                               srcdir)
        gencpy(srcdir,destdir,sh_pp_fyle)

        if not os.path.exists(srcdir+'/' + sh_md_fyle):
            raise RuntimeError(sh_md_fyle+" not found in "+\
                               srcdir)
        gencpy(srcdir,destdir,sh_md_fyle)
        

    if continue_run == 1: #at least one tpr file present
        sh_mdrun = sh_md_fyle.replace('_pyinp','')
        print("Recopying ", sh_md_fyle)
        gencpy(srcdir,destdir,sh_md_fyle)

    return continue_run
#------------------------------------------------------------------
# Edit pre-processing shell script files
def edit_sh_files(jtype,inp_type,polycfg,topfyle,sh_fyle,tempval,\
                  maindir,destdir):

    if not os.path.exists(maindir + '/' + sh_fyle):
        raise RuntimeError(sh_fyle + ' not found in ' + maindir)

    dval = 0.3
    
    # job/box name
    jname = jtype+'_'+inp_type+ '_T' + str(tempval)
    box_conffyle = "initconf.gro"

    # edit sh_fyle
    py_fname = sh_fyle
    rev_fname = py_fname.replace('_pyinp','')
    fr  = open(maindir + '/' + py_fname,'r')
    fw  = open(destdir + '/' + rev_fname,'w')
    fid = fr.read().replace("py_jobname",jname).\
          replace("py_meltconf",polycfg).\
          replace("py_boxmeltconf",box_conffyle).\
          replace("py_topol",topfyle).\
          replace("py_finconf",box_conffyle).\
          replace("py_dval",str(dval))
    fw.write(fid)
    fr.close(); fw.close()
#------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------

