# Generate initial structures/patches for modified cellulose systems
# Reference: Loukas's code for cellulose build
# Version: Spet-14-2021

# Generates single monomer with acetylated patch for checking

# Load packages
package require psfgen

# Set inputs
set deg_poly 1 ; # DONT CHANGE THIS
set carb_pos 6 ; #carbon position for substitution

# Define variables
set outname singlemon_6AC5
set dirname ../elementary_fibrils

# Load topology
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_carb.rtf
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/toppar_all36_carb_imlab.str
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all35_ethers.rtf

# Build acetylated monomer
set chain C 
foreach n  { 0 } {
    resetpsf    
    set segment $chain$n
    segment  $segment {
	
	for {set i 1} {$i<=$deg_poly} {incr i} {
	    residue $i BGLC 
	}
    }
    patch 6AC5 $segment:1
    
    coordpdb ../elementary_fibrils/csff-elementary-fibril-$segment.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf
    unset segment
    resetpsf
}



#------------Combine psf/pdb---------------------------------------
resetpsf
foreach i  { 0 } {
    readpsf  C$i.psf
    coordpdb C$i.pdb
}

guesscoord
writepdb $outname.pdb	 
writepsf $outname.psf

exit
