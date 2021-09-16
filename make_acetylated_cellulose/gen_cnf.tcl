# Generate initial structures/patches for modified cellulose systems
# Reference: Loukas's code for cellulose build
# Version: Feb-18-2021

# Note: puts command to fp will write to a separate file for
# checking. Please comment out if not needed.

# Load packages
package require psfgen

# Set inputs
set deg_poly 5
set num_chains 5
set acet_prob 0.3

# Define variables
set chnamepref ace ; # chain name prefix
set outname modifiedCNF_m1
set inpname acetyl_cnf_m1

# Load topology
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_carb.rtf
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_cgenff.rtf
topology /home/v0e/allcodes/files_cnf/acetylation/m1/$inpname.rtf

# Set random number generator
set t1 [clock milliseconds]
expr { srand(${t1}) }

# Open structure file for cross checking
set fp [open "$outname.tcl" w]

# Begin building chains
for {set chcntr 1} {$chcntr <= $num_chains} {incr chcntr 1} {
    puts $fp "# Chain number: $chcntr"
    puts $fp "resetpsf"
    set segment $chnamepref$chcntr
    resetpsf
    puts $fp "segment $segment {"
    segment $segment {
	#set first monomer to be BGLC since that is the guess input
	puts $fp "residue 1 BGLC" 
	residue 1 BGLC
	for {set i 2} {$i <= $deg_poly} {incr i 1} {
	    if {rand() < $acet_prob} {
		puts $fp "residue $i BGLC"
		residue $i BGLC
	    } else {
		puts $fp "residue $i lig"
		residue $i LIG
	    }
	}
    }
    puts $fp "}"
    puts $fp "# Set patches"
    for {set i $deg_poly} {$i > 1} {incr i -1} {
	set j [expr $i-1]
	patch 14bb $segment:$i $segment:$j
	puts $fp "patch 14bb $segment:$i $segment:$j"
    }
    coordpdb C0.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf
    unset segment
    resetpsf
}

# Combine psf/pdb
resetpsf
for {set chcntr 1} {$chcntr <= $num_chains} {incr chcntr 1} {
    readpsf  $chnamepref$chcntr.psf
    coordpdb $chnamepref$chcntr.pdb
}

guesscoord
writepdb $outname.pdb	 
writepsf $outname.psf

exit
