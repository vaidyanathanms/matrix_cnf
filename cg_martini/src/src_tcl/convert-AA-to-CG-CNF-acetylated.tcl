### map the atomistic structure to martini beads for acetylated cellulose

package require mol
set nchains 18
set DP 20
# for any given monomer without acetylation, the maximum atoms per resid cannot exceed 23, when it does check for acetylated atom types 
set max_atoms_per_res 23

set m_length [lindex $argv 0]
### get all the files for a given acetylated composition
set dir "."
set files [glob -directory $dir initstruct*.pdb]
set files [lsort -dictionary $files]
set counter 1

foreach file $files {
  mol load pdb $file
# Set bead types for acetylated segment (B4)
  for {set frag 0} {$frag < $nchains} {incr frag 1} {
   for {set res 1} {$res <= $DP} {incr res 1} {
   set indice_lst [atomselect top "fragment $frag and resid $res"]
   set indices [eval $indice_lst num]
   if { $indices > $max_atoms_per_res} {
     if {  $m_length eq 1 } {
        set bx [atomselect top "(name CY6 CX6 OX6) and fragment $frag and resid $res"]
        puts "[eval $bx num]"
        set comx [measure center $bx weight mass]
        set cx [atomselect top "name CY6 and fragment $frag and resid $res"]
        $cx moveto $comx
        $cx set name B4
      } elseif { $m_length eq 3 } {
        set bx [atomselect top "(name CV6 CX6 OX6) and fragment $frag and resid $res"]
        set comx [measure center $bx weight mass]
        set by [atomselect top "(name CW6 CY6) and fragment $frag and resid $res"]
        set comy [measure center $by weight mass]
        set cx [atomselect top "name CY6 and fragment $frag and resid $res"]
        $cx moveto $comx
        $cx set name B4
        set cy [atomselect top "name CV6 and fragment $frag and resid $res"]
        $cy moveto $comy
        $cy set name B5      
      }
     }
   }
  }
  
  # Set bead types B1 and B3 after finding COM (P4 and P1)
  
  for {set frag 0} {$frag < $nchains} {incr frag 1} {
   for {set res 1} {$res <= $DP} {incr res 1} {
   set b3 [atomselect top "(name O6 C6 C5) and fragment $frag and resid $res"]
   set com1 [measure center $b3 weight mass]
   set c6 [atomselect top "name C6 and fragment $frag and resid $res"]
   $c6 moveto $com1
   $c6 set name B3
   set b1 [atomselect top "(name C2 C3 O2 O3) and fragment $frag and resid $res"]
   set com2 [measure center $b1 weight mass]
   set o2 [atomselect top "name O2 and fragment $frag and resid $res"]
   $o2 moveto $com2
   $o2 set name B1
   }
  }
  # Set bead type B2 for all residues except 1 and 20 in each chain after finding COM (PX)
  
  for {set frag 0} {$frag < $nchains} {incr frag 1} {
   for {set res 2} {$res <= $DP} {incr res 1} {
   set b2 [atomselect top "(name C4 O5 C1 O4) and fragment $frag and resid $res"]
   set com1 [measure center $b2 weight mass]
   set c4 [atomselect top "name C4 and fragment $frag and resid $res"]
   $c4 moveto $com1
   $c4 set name B2
   }
  }
  
  # Set bead type B2 for residue 1 in each chain after finding COM (PX)
  
  for {set frag 0} {$frag < $nchains} {incr frag 1} {
   for {set res 1} {$res < 2} {incr res 1} {
   set b2 [atomselect top "(name O4 C4 O5 C1 O4) and fragment $frag and resid $res"]
   set com1 [measure center $b2 weight mass]
   set c4 [atomselect top "name C4 and fragment $frag and resid $res"]
   $c4 moveto $com1
   $c4 set name B2
   }
  }
  
  # Set bead type B2 for residue 20 in each chain after finding COM (PX)
  
  for {set frag 0} {$frag < $nchains} {incr frag 1} {
   for {set res $DP} {$res <= $DP} {incr res 1} {
   set b2 [atomselect top "(name C4 O5 C1 O1) and fragment $frag and resid $res"]
   set com1 [measure center $b2 weight mass]
   set c4 [atomselect top "name C4 and fragment $frag and resid $res"]
   $c4 moveto $com1
   $c4 set name B2
   }
  }  
 
  if {  $m_length eq 1 } { 
     set towrite [atomselect top "name B1 B2 B3 B4"]
  } elseif {  $m_length eq 3 } {
     set towrite [atomselect top "name B1 B2 B3 B5 B4"]
  }
     $towrite writepdb test.pdb
  
  mol load pdb test.pdb
  
  set existing_B3type [atomselect top "name B3"]
  $existing_B3type set name P5 
  	
  set existing_B1type [atomselect top "name B1"]
  $existing_B1type set name P6
  
  set existing_B2type [atomselect top "name B2"]
  $existing_B2type set name PX

  if {  $m_length eq 1 } {
   set existing_B4type [atomselect top "name B4"]
   $existing_B4type set name SN4a 
  } elseif {  $m_length eq 3 } {
  set existing_B4type [atomselect top "name B4"]
  $existing_B4type set name SN4a
  set existing_B5type [atomselect top "name B5"]
  $existing_B5type set name TC1
    
  }  
  
  set sel_all [atomselect top all]
  $sel_all set resname CEL1X
  set new_file [ string trim $file "-AA.pdb"]
  puts "######################*****"
  $sel_all writepdb .$new_file-Martini.pdb
  $sel_all writegro .$new_file-Martini.gro
  incr counter 1
  rm test.pdb
}

exit
