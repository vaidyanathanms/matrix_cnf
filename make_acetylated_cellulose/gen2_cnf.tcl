# Generate initial structures/patches for modified cellulose systems
# Reference: Loukas's code for cellulose build
# Version: June-01-2021

# Uses 6TAC patch from /toppar_all36_carb_imlab.str and BGLC from
# top_all36_carb.rtf 
# This method was easier than creating one acetylated carbon using
# ligand builder and matching the carbons.

# Note: puts command to fp will write to a separate file for
# checking structures. Please comment out if not needed.

# Load packages
package require psfgen

# Set inputs
set deg_poly 50
set num_chains 5
set acet_prob 0.3
set carb_pos 6 ; #carbon position for substitution

# Since first monomer is cellulose, the acetylation probability needs
# to be renormalized
set new_acet_prob [ expr ($acet_prob*$deg_poly)/($deg_poly-1) ]

# Define variables
set chname1 AC$carb_pos ; # chain name prefix for chain ID <= 9
set chname2 A$carb_pos  ; # chain name prefix for chain ID > 9
set outname modifiedCNF_m1
set inpname acetyl_cnf_m1
set dirname ../elementary_fibrils

# Load topology
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_carb.rtf
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/toppar_all36_carb_imlab.str

# Set random number generator
set t1 [clock milliseconds]
expr { srand(${t1}) }

# Open structure file for cross checking
set fp [open "$outname.tcl" w]
puts $fp "Old acetylation probability: $acet_prob"
puts $fp "New acetylation probability: $new_acet_prob"

# Begin building center chains
for {set chcntr 1} {$chcntr <= $num_chains} {incr chcntr 1} {

    puts $fp "# Chain number: $chcntr"
    # Maximum length of segment number is 4
    if {$chcntr > 10} {
	set segment $chname1$chcntr
    } else {
	set segment $chname2$chcntr
    }

    puts $fp "resetpsf" 
    resetpsf
    puts $fp "segment $segment {"
    segment $segment {
	# set first monomer as BGLC since it is obtained from C0.pdb
	# of cellulosic structure. Assign random numbers during
	# patches to decide whether the monomer is cellulose or
	# acetylated-cellulose. 

	for {set i 1} {$i <= $deg_poly} {incr i 1} {
	    residue $i BGLC
	    puts $fp "residue $i BGLC" 
	}
    }
    puts $fp "}"

    puts $fp "# Set patches"
    for {set i $deg_poly} {$i > 1} {incr i -1} {
	set j [expr $i-1]
	# Convert to acetylated-cellulose if the random number is less
	# than the acet_prob. DO NOT conver the first monomer. No need
	# of an if statement since by decrementing, we are looking at
	# values of i>1.
	if {rand() < $new_acet_prob} {
	    patch 6TAC $segment:$i 
	    puts $fp "patch 6TAC $segment:$i"
	} 
	patch 14bb $segment:$i $segment:$j
	puts $fp "patch 14BB $segment:$i $segment:$j"
    }
    coordpdb ./start_struct/C0.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf

    # Unset everything
    unset segment
    resetpsf
}

# Combine psf/pdb
resetpsf
for {set chcntr 1} {$chcntr <= $num_chains} {incr chcntr 1} {
    if {$chcntr > 10} {
	readpsf  $chname1$chcntr.psf
	coordpdb $chname1$chcntr.pdb
    } else {
	readpsf  $chname2$chcntr.psf
	coordpdb $chname2$chcntr.pdb
    }	
}

guesscoord
writepdb $outname.pdb	 
writepsf $outname.psf

exit
