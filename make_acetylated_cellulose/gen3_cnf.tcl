# Generate initial structures/patches for modified cellulose systems
# Reference: Loukas's code for cellulose build
# Version: June-10-2021

# Uses 6TAC patch from /toppar_all36_carb_imlab.str and BGLC from
# top_all36_carb.rtf 
# This method was easier than creating one acetylated carbon using
# ligand builder and matching the carbons.

# Note: puts command to fp will write to a separate file for
# checking structures. Please comment out if not needed.

# Main change from Version 2: Use center and origin chains from
# Loukas's cellulose code and then patch with 6TAC. This way there is
# no need to renormalize the probability



# Load packages
package require psfgen

# Set inputs
set deg_poly 20
set acet_prob 0.3
set carb_pos 6 ; #carbon position for substitution

# Define variables
set outname modifiedCNF_m1
set dirname ../elementary_fibrils
set acet_cnt 0

# Load topology
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_carb.rtf
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/toppar_all36_carb_imlab.str

# Set random number generator
set t1 [clock milliseconds]
expr { srand(${t1}) }

# Open structure file for cross checking
set fp [open "acetylation.dat" w]
puts $fp "Acetylation probability: $acet_prob"

#----------------Non-acetyated chains----------------------------------
# Build core center chains (no acetylation)
set chain C
foreach n { 1 4 5 6 10 } {
    resetpsf
    set segment $chain$n
    segment  $segment {
	for {set i 1} {$i<=$deg_poly} {incr i} {
	    residue $i BGLC 
	}
    }
    for {set i $deg_poly} {$i>1} {incr i -1} {
	set j [ expr $i -1]
	patch   14bb $segment:$i $segment:$j
    }
 
    coordpdb ../elementary_fibrils/csff-elementary-fibril-$segment.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf
    unset segment
    resetpsf
}

# Build core origin chains (no acetylation)
set chain O
foreach n { 1 2 6 7 } {
    resetpsf
    set segment $chain$n
    segment  $segment {
	for {set i 1} {$i<=$deg_poly} {incr i} {
	    residue $i BGLC 
	}
    }
    for {set i $deg_poly} {$i>1} {incr i -1} {
	set j [ expr $i -1]
	patch   14bb $segment:$i $segment:$j
    }
 
    coordpdb ../elementary_fibrils/csff-elementary-fibril-$segment.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf
    unset segment
    resetpsf
}

#-Acetylated chains: acetlyate equally among different chains------

# Build acetylated center chains

set chain C 
foreach n  { 0 2 7 9 11 } {
    resetpsf    
    set segment $chain$n
    puts $fp "# Chain ID: $segment"
    segment  $segment {
	for {set i 1} {$i<=$deg_poly} {incr i} {
	    residue $i BGLC 
	}
    }
    for {set i $deg_poly} {$i>1} {incr i -1} {
	set j [ expr $i -1]
	patch   14bb $segment:$i $segment:$j
	if {rand() < $acet_prob} {
	    patch 6TAC $segment:$i 
	    puts $fp "patch 6TAC $segment:$i"
	    set acet_cnt [ expr $acet_cnt + 1 ]
	} 
    }
    coordpdb ../elementary_fibrils/csff-elementary-fibril-$segment.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf
    unset segment
    resetpsf
}


# Build acetylated origin chains
set chain O
foreach n { 0 3 5 8 } {
    resetpsf
    set segment $chain$n
    puts $fp "# Chain ID: $segment"
    segment  $segment {
	for {set i 1} {$i<=$deg_poly} {incr i} {
	    residue $i BGLC 
	}
    }
    for {set i $deg_poly} {$i>1} {incr i -1} {
	set j [ expr $i -1]
	patch   14bb $segment:$i $segment:$j
	if {rand() < $acet_prob} {
	    patch 6TAC $segment:$i 
	    puts $fp "patch 6TAC $segment:$i"
	    set acet_cnt [ expr $acet_cnt + 1 ]
	}
    }
 
    coordpdb ../elementary_fibrils/csff-elementary-fibril-$segment.pdb $segment
    guesscoord
    regenerate angles dihedrals
    writepdb $segment.pdb
    writepsf $segment.psf
    unset segment
    resetpsf
}
puts $fp "Total number of acetylated monomers: $acet_cnt"
set sim_prob [ expr ($acet_cnt)/(18.0*$deg_poly) ]
puts $fp "Simulated acetylation: $sim_prob"
#------------Combine psf/pdb---------------------------------------
resetpsf
foreach i  { 1 2 6 7 } {
    readpsf  O$i.psf
    coordpdb O$i.pdb
}
foreach i  { 1 4 5 6 7 10 } {
    readpsf  C$i.psf
    coordpdb C$i.pdb
}
foreach i  { 0 3 5 8 } {
    readpsf  O$i.psf
    coordpdb O$i.pdb
}
foreach i  { 0 2 9 11 } {
    readpsf  C$i.psf
    coordpdb C$i.pdb
}

guesscoord
writepdb $outname.pdb	 
writepsf $outname.psf

exit
