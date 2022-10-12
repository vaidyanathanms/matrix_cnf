# User needs to input data for all this case
# Input/Output directory paths

acet_dir   = '/home/v0e/allcodes/files_cnf/main/acetylation_files' # input dir for acetylation
cell_topdir= '/home/v0e/allcodes/files_cnf/main/cell_toppar' # top/par dir for cellulose/acetylated cell
natv_cnfdir= '/home/v0e/allcodes/files_cnf/main/elementary_fibrils' #cellulose inps
poly_mat   = '/home/v0e/allcodes/files_cnf/polymer_matrices' #i/o dir poly matrices
pack_exe   = '/home/v0e/packmol/packmol' # packmol executable
mdp_dir    = '/home/v0e/allcodes/files_cnf/main/mdp_files' # mdp files for GMX
scr_dir    = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir (MD run dir)
#------------------------------------------------------------------
# NOTE: Directory from CHARMM should be inside "poly_mat" (see above)# directory with the prefix charmm_
# For instance, if PLA is the matrix, then the files from CHARMM
# should be in a folder named charmm_pla and the folder should be in
# poly_mat/charmm_pla, where poly_mat is the variable above.
