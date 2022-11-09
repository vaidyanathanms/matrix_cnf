# This is to generate a topology file for GROMACS
# Uses writegmxtop from topotools
# Note 1: This uses a dummy file of the modified CNF and 
# Note 2: Uses the modified CNF pdb/psf files
# Note 3: Make sure no unnecessary bonds are generated

# Load packages
package require psfgen
package require topotools

# Load topology files
topology ../cell_toppar/top_all36_carb.rtf ;# for basic cellulose 
topology ../cell_toppar/toppar_all36_carb_imlab.str ;# for modifications

# Read modified CNF pdb/psf files
resetpsf
mol new  modified_m1.psf
mol addfile modified_m1.pdb

# Generate GROMACS *.top file
topo writegmxtop modified_m1.top [list ../cell_toppar/par_all36_carb.prm ../cell_toppar/toppar_all36_carb_imlab.str]

exit



