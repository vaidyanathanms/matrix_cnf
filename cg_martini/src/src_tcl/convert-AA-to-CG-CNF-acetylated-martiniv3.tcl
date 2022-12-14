### map the atomistic structure to martini beads for acetylated cellulose

package require mol
######## Read the input coordinate file in pdb/gro format
set file_type [lindex $argv 0]
set file_string [lindex $argv 1]
puts $file_string
set nchains [lindex $argv 2] 
set DP [lindex $argv 3]
set cnfs_per_bundle [lindex $argv 4]
set m_length [lindex $argv 5]
set nfragments [expr $nchains * $cnfs_per_bundle]
set nsegments [expr $nchains * $cnfs_per_bundle]

# for any given monomer without acetylation, the maximum atoms per residue cannot exceed 23, when it does check for acetylated atom types 
set max_atoms_per_res 23

### get all the files for a given acetylated composition
set dir "."
set coord_files [glob -directory $dir initstruct*]
set coord_files [lsort -dictionary $coord_files]
set counter 1

foreach file $coord_files {
  if { $file_type eq "pdb" } {
        mol load pdb $file
     } elseif { $file_type eq "gro" } {
        mol load gro $file
     } else {
        puts "incorrect file type"
        exit
     }

# Set bead types for acetylated segment (B5 and/or B6)
  for {set frag 0} {$frag < $nfragments} {incr frag 1} {
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
        $cx set name B5
      } elseif { $m_length eq 3 } {
        set bx [atomselect top "(name CV6 CX6 OX6) and fragment $frag and resid $res"]
        set comx [measure center $bx weight mass]
        set by [atomselect top "(name CW6 CY6) and fragment $frag and resid $res"]
        set comy [measure center $by weight mass]
        set cx [atomselect top "name CY6 and fragment $frag and resid $res"]
        $cx moveto $comx
        $cx set name B5
        set cy [atomselect top "name CV6 and fragment $frag and resid $res"]
        $cy moveto $comy
        $cy set name B6      
      }
     }
   }
  }
  for {set frag 0} {$frag < $nfragments} {incr frag 1} {
   for {set res 1} {$res <= $DP} {incr res 1} {
      set b1 [atomselect top "(name O6 HO6 C6 H61 H62) and fragment $frag and resid $res"]
      set com1 [measure center $b1 weight mass]
      if {$res eq 1} {
         set b2 [atomselect top "(name O5 C1 H1 O1 HO1) and fragment $frag and resid $res"]
         set com2 [measure center $b2 weight mass]
       } else {
         set b2 [atomselect top "(name O5 C1 H1) and fragment $frag and resid $res"]
         set com2 [measure center $b2 weight mass]
       }
      set b3 [atomselect top "(name C2 H2 C3 H3 O2 HO2 O3 HO3) and fragment $frag and resid $res"]
      set com3 [measure center $b3 weight mass]
      if {$res eq $DP} {
         set b4 [atomselect top "(name C5 H5 C4 H4 O4 HO4) and fragment $frag and resid $res"]
         set com4 [measure center $b4 weight mass]
       } else {
         set b4 [atomselect top "(name C5 H5 C4 H4 O4) and fragment $frag and resid $res"]
         set com4 [measure center $b4 weight mass]
       }
      set c1 [atomselect top "name C1 and fragment $frag and resid $res"]
      $c1 moveto $com1
      $c1 set name B1
      set c2 [atomselect top "name C2 and fragment $frag and resid $res"]
      $c2 moveto $com2
      $c2 set name B2
      set c3 [atomselect top "name C3 and fragment $frag and resid $res"]
      $c3 moveto $com3
      $c3 set name B3
      set c4 [atomselect top "name C4 and fragment $frag and resid $res"]
      $c4 moveto $com4
      $c4 set name B4
   }
  } 
  if {  $m_length eq 1 } { 
     set towrite [atomselect top "name B1 B2 B3 B4 B5"]
  } elseif {  $m_length eq 3 } {
     set towrite [atomselect top "name B1 B2 B3 B4 B5 B6"]
  }
  $towrite writepdb test.pdb
  mol load pdb test.pdb
  set T1 [atomselect top "name B1"]
  $T1 set name T1
  set T2 [atomselect top "name B2"]
  $T2 set name T2
  set R3 [atomselect top "name B3"]
  $R3 set name R3
  set S4 [atomselect top "name B4"]
  $S4 set name S4
  puts "print"
  if {  $m_length eq 1 } {
   set A1 [atomselect top "name B5"]
   $A1 set name A1
  } elseif {  $m_length eq 3 } {
  set A1 [atomselect top "name B5"]
  $A1 set name A1
  set A2 [atomselect top "name B6"]
  $A2 set name A2  
  }  
  set sel [atomselect top "resid 20 and name T2"]
  $sel set name S2
  set sel_all [atomselect top all]
  $sel_all set resname CELL
  set new_file [ string trim $file "-AA.pdb"]
  puts "######################*****"
  $sel_all writepdb .$new_file-Martini.pdb
  $sel_all writegro .$new_file-Martini.gro
  incr counter 1
  rm test.pdb
}

exit
