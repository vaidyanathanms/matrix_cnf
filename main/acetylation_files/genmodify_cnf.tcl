# Generate initial structures/patches for modified cellulose systems
# Generalized python code
# Reference: Loukas's code for cellulose build
# Version: Sept-01-2021

# Uses acetyle patch from /toppar_all36_carb_imlab.str and BGLC from
# top_all36_carb.rtf 
# Acetyl patches: TAC6 (m1), 6AC11 (m3), 6AC11 (m7), 6A11 (m11)
# This method was easier than creating one acetylated carbon using
# ligand builder and matching the carbons.

# Note: puts command to fp will write to a separate file for
# checking structures. Please comment out if not needed.

# Main change from Version 2: Use center and origin chains from
# Loukas's cellulose code and then use appropriate acetyl patch. 
# This way there is no need to renormalize the probability.

# Difference from V3.0: acetylates only exposed O6

# Load packages
package require psfgen

# Set inputs
set deg_poly 20
set acet_prob 0.3
set carb_pos 6 ; #carbon position for substitution
set targ_relerr 0.1 ; # value between 0 and 1
set max_att 100 ; # maximum attempt for acetylation

# Define variables
set outname /home/v0e/allcodes/files_cnf/modified_m11
set dirname ../elementary_fibrils

# Load topology
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_carb.rtf
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/toppar_all36_carb_imlab.str

# Set random number generator
set t1 [clock milliseconds]
expr { srand(${t1}) }

# Open structure file for cross checking
set fp [open "acetylation.dat" w]
puts $fp "Acetylation probability: $acet_prob"
puts $fp "Acetylation patch used: 6AC11"

#----------------Non-acetyated chains----------------------------------
# Build core center chains (no acetylation)
set chain C
foreach n { 4 5 6} {
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
set trialnum 0
set rel_err 100
set max_acet [ expr ((11.0*$deg_poly)/(2.0)) ]
puts $fp "Possible number of monomers acetylated: $max_acet"
set sim_poss [ expr double((round($acet_prob*$max_acet))) ]
puts $fp "Weighted possible number of acetylated monomers: $sim_poss"

while { ($rel_err > $targ_relerr) || ($trialnum > $max_att) } {

    set trialnum [ expr ($trialnum + 1) ]
    puts "Current relative error: $rel_err"
    puts "Trial number: $trialnum"
    puts $fp "********Trial: $trialnum*******************************"
    set acet_cnt 0

    # Build acetylated center chains
    set chain C 
    foreach n  { 0 1 2 7 9 10 11 } {
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
	    # odd monomers to be acetylated
	    if {$n == 0 | $n == 9 | $n == 1} {
		if {$i%2 == 1} { 
		    if {rand() < $acet_prob} {
			patch 6AC11 $segment:$i 
			puts $fp "patch 6AC11 $segment:$i"
			set acet_cnt [ expr $acet_cnt + 1 ]
		    } 
		}
	    } elseif {$n == 2 | $n == 7 | $n == 11 | $n == 10} {
		if {$i%2 == 0} { 
		    if {rand() < $acet_prob} {
			patch 6AC11 $segment:$i 
			puts $fp "patch 6AC11 $segment:$i"
			set acet_cnt [ expr $acet_cnt + 1 ]
		    } 
		}
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
	    # odd monomers to be acetylated
	    if {$n == 0 | $n == 5} {
		if {$i%2 == 1} { 
		    if {rand() < $acet_prob} {
			patch 6AC11 $segment:$i 
			puts $fp "patch 6AC11 $segment:$i"
			set acet_cnt [ expr $acet_cnt + 1 ]
		    } 
		}
	    } elseif {$n == 3 | $n == 8} {
		if {$i%2 == 0} { 
		    if {rand() < $acet_prob} {
			patch 6AC11 $segment:$i 
			puts $fp "patch 6AC11 $segment:$i"
			set acet_cnt [ expr $acet_cnt + 1 ]
		    } 
		}
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
set rel_err  [ expr abs( (double($acet_cnt) - double($sim_poss))/(double($sim_poss)) ) ]
puts $fp "Relative error: $rel_err"
puts $fp "**********************************************************"
}

if {$rel_err < $targ_relerr} {
    puts $fp "Final configuration obtained"
} else {
    puts $fp "ERR: No final configuration obtained after $max_att trials"
}

#------------Combine psf/pdb---------------------------------------
resetpsf
foreach i  { 1 2 6 7 } {
    readpsf  O$i.psf
    coordpdb O$i.pdb
}
foreach i  { 4 5 6 } {
    readpsf  C$i.psf
    coordpdb C$i.pdb
}
foreach i  { 0 3 5 8 } {
    readpsf  O$i.psf
    coordpdb O$i.pdb
}
foreach i  { 0 1 2 7 9 10 11 } {
    readpsf  C$i.psf
    coordpdb C$i.pdb
}

guesscoord
writepdb $outname.pdb	 
writepsf $outname.psf

exit
