# To check guassian polymer matrix chains
# Add above to PACKMOL along with modified cellulose chains or additional
# polymers

import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
import dir_names
import mdp
from aux_pack import * # function definitions
#------------------------------------------------------------------

# Polymer matrix
matrix   = 'pla' #pla/pp/petg/p3hb

# Input data - Polymer matrix
mat_pdb  = 'step3_input.pdb' # matrix input pdb file - ONLY PDB
nmons    = 40 # number of matrix monomers in output
nchains  = 81 # number of matrix chains in output
gaus_tol = 0.05 # tolerance for checking gaussianity

# Input data - Cellulose/Acetylated Cellulose
acet_val = 1 # m1 - 1, m3 - 3, m7 - 7, m11 - 11
acet_per = 0.5 # fraction of acetylated cellulose
acet_new = 1 # 0 - use old, 1-delete and regenerate
ncnf     = 1 # number of cellulose bundles (18 chains)
cell_dp  = 20 # degree of polymerization of cellulose chains
acet_tol = 0.1 # tolerance for acetylation
acetpref = 'modified_m' # prefix for acetylated files
acet_att = 100 # maximum attempts to create acetylated cellulose

# Input data - Additives - blends/triblock copolymers
add_poly = 'None' # blend/triblock
ex_ptype = ['paa','pvp'] #p1,p2 for blend, p1_p2_p3 for triblock
ex_nch   = [9, 9] # number of chains of each type
ex_nmon  = [18, 11] # degree of polymerization of each type

# Input data - Packmol
inppack  = 'pack_cellulose.inp' # PACKMOL input file
fin_box  = 1.1 # final box size relative to max dimension of cnf/matrix
run_pack = 1 # 1-run packmol
cleandir = 1 # clean directories
packsh   = 'run_packmol_pyinp.sh'

# Input data - GROMACS
set_mdp   = 1 # Copy mdp files
mdp_files = ['minim_pyinp.mdp','nvt_pyinp.mdp','nvt_high_pyinp.mdp'\
             ,'npt_berendsen_pyinp.mdp','npt_main_pyinp.mdp']
Thigh  = 600 # High temperature for equilibration
Ttarg  = 300 # Target temperature for simulations
refP   = 1   # Reference pressure
pp_run = 'run_preprocess_pyinp.sh'
md_run = 'run_md_pyinp.sh'
#------------------------------------------------------------------

# Import directory paths
main_dir    = os.getcwd() # current dir
acet_dir    = dir_names.acet_dir
cell_topdir = dir_names.cell_topdir
natv_cnfdir = dir_names.natv_cnfdir
poly_mat    = dir_names.poly_mat
pack_exe    = dir_names.pack_exe
mdp_dir     = dir_names.mdp_dir
scr_dir     = dir_names.scr_dir

# Check for directory paths and input consistency
chrm_dir = poly_mat + '/charmm_' + matrix #CHARMM inp dir for poly matrices
pack_sup = poly_mat + '/cnf_packed_' + matrix + '_N_' + str(nmons)\
           + '_M_' + str(nchains)  # final packed output super dir
gmx_mat  = chrm_dir + '/gromacs'

check_dir(natv_cnfdir)
check_dir(cell_topdir)
check_dir(poly_mat)
check_dir(chrm_dir)
check_dir(gmx_mat)
check_dir(scr_dir)

if fin_box < 1:
    raise RuntimeError("Unphysical final box magnification\n")
#------------------------------------------------------------------

# Make output directories/files
scr_dir   = scr_dir + '/cellulose_nanofibers'
if acet_per != 0:
    acet_fyle = acetpref + str(acet_val)
    mod_cell  = 1
    acet_outdir = '../../acetyl_cellulose' #out acet_cell dir
else:
    acet_fyle = native_cell_dp + str(cell_dp)
    mod_cell = 0
#------------------------------------------------------------------

# Main code

print('Begin creating packmol files...')
print('Begin checking input files and create output directories...')

# Check for polymer matrix files 
check_inp_files(gmx_mat,mat_pdb) 

# Create output directories
pack_dir = create_output_dirs(pack_sup,acet_val,acet_per,add_poly)

print('Checking matrix input files for gaussian chains...')
# Check for Gaussian input chains
polygausdir,rgmax = check_gaussianity_and_write(gmx_mat,mat_pdb,nmons,\
                                                nchains,matrix,gaus_tol,pack_dir)

# Acetylate chains if needed
if acet_per != 0:
    print('Making acetylated chains ..')
    make_acet_cell(acet_dir,acet_val,cell_dp,acet_per,acet_tol,\
                   acet_fyle,acet_att,pack_dir,acet_new,natv_cnfdir)
    os.chdir(main_dir)

print('Determining optimal dimensions of the box...')
# Check for cnf dimensions and compare with rg of polymer matrix
xmin,ymin,zmin,xmax,ymax,zmax = find_cnf_dim(acet_dir,acet_fyle,pack_dir)
dmax = max(xmax-xmin,ymax-ymin,zmax-zmin,rgmax)

# Start writing PACKMOL files
packfyle = pack_dir + '/' + inppack
fpack   = open(packfyle,'w')
packed_cnfpdb = packmol_headers(fpack,matrix,pack_dir,\
                                acet_per,acet_val,add_poly)

# Pack CNF/matrix/additional polymers
print('Writing packmol scripts for packing cellulose and matrix chains...')
pack_cellulose_chains(fpack,pack_dir,acet_fyle,ncnf,xmin,ymin,zmin,xmax,ymax,zmax,fin_box,dmax)
pack_polymer_matrix(polygausdir,matrix,nchains,xmin,ymin,zmin,xmax,ymax,zmax,fpack,fin_box,dmax)
if add_poly.lower() != 'None'.lower():
    exgausdir,exrgmax = check_gaussianity_and_write(gmx_mat,mat_pdb,nmons,\
                                                    nchains,matrix,gaus_tol,\
                                                    pack_dir)
    pack_extra_poly(exgausdir,matrix,ex_nch,xmin,ymin,zmin,xmax,ymax,zmax,fpack,fin_box,dmax)

fpack.close() # close PACKMOL input file

# Run PACKMOL
if run_pack:
    print('Submitting PACKMOL scripts..')
    run_packmol(packfyle,pack_exe,pack_dir,packsh,main_dir)

# top/prm/itp file arrays
prmfyle_arr = []; itpfyle_arr = []

# Generate and split top file for cell/acetylated cellu
print('Generating GMX top files for celluloses..')
make_top_file_for_acetcell(main_dir,pack_dir,cell_topdir,acet_fyle)
prm_file,itp_file,mol_infoarr = split_top_file_to_prmitp(acet_fyle,pack_dir,main_dir)
prmfyle_arr.append(prm_file); itpfyle_arr.append(itp_file)

# Copy itp/top/prm file for polymer matrix
# Check for toppar inside gromacs directory output by CHARMM-GUI
print('Copying toppar files from CHARMM-GUI for polymer matrices.')
copy_mat_toppar_files(gmx_mat,pack_dir,prmfyle_arr,\
                      itpfyle_arr,mol_infoarr)

# Add ; to all the polymer matrix forcefield files (prm)
print('Editing forcefield files of polymer matrices..')
add_comment_to_ff_files(prmfyle_arr)

# Combine polymer matrix and acet cell files into one top file
print('Combining prm/itp files of cellulose and polymers'\
      + ' into one single top file for GMX calculations..')
out_topo_file = combine_top_files(pack_dir,prmfyle_arr,\
                                  itpfyle_arr,mol_infoarr,\
                                  nchains)

# Clean up psf/pdb files of native cellulose
clean_and_sort_files(pack_dir,acetpref)

# Create final directories in scratch
if set_mdp:
    check_dir(mdp_dir)
    print('Copying  directories of MD simulations...')
    # Set thermostat/top variables
    Tetau_nvt,Tetau_highnvt,Tetau_berend,Tetau_parrah,\
        Prtau_berend,Prtau_parrah,ref_pres,\
        melt_topname=mdp.couple_coeff('melts','None')

    # Find tc_groups
    # Copy and edit mdp files (all temperatures are same)
    mdp.check_cpy_mdp_files(mdp_dir,pack_dir,mdp_files,'melts'\
                            ,Tetau_nvt,Tetau_highnvt,Tetau_berend\
                            ,Tetau_parrah,Prtau_berend,Prtau_parrah\
                            ,Ttarg,Thigh,refP,'System','Single',\
                            main_dir,'None')


    # Generate shell input files
    mdp.edit_sh_files('pp',matrix,packed_cnfpdb,out_topo_file,\
                      pp_run,Ttarg,mdp_dir,pack_dir)#preprocess
    mdp.edit_sh_files('md',matrix,packed_cnfpdb,out_topo_file,\
                      md_run,Ttarg,mdp_dir,pack_dir)#md run

os.chdir(main_dir)
    
print('Done :) ..')



