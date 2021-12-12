# This is to generate a topology file for cellulose/acetylated cellulose in GMX format
# Uses writegmxtop from topotools
# Note 1: This uses a dummy file of the modified CNF and 
# Note 2: Uses the (un)modified CNF pdb/psf files
# Note 3: Make sure no unnecessary bonds are generated

# Load packages
package require psfgen
package require topotools

# Load topology files
topology py_topdir/top_all36_carb.rtf ;# for basic cellulose 
topology py_topdir/toppar_all36_carb_imlab.str ;# for modifications

# Read (un)modified CNF pdb/psf files
resetpsf
set cellfile py_outname
mol new  $cellfile.psf
mol addfile $cellfile.pdb

# Generate GROMACS *.top file
topo writegmxtop $cellfile.top [list py_topdir/par_all36_carb.prm py_topdir/toppar_all36_carb_imlab.str]

exit



