package require mol

set dir [lindex $argv 0]
puts "$dir"
set files [glob -directory "$dir" *.pdb]
set files [lsort -dictionary $files]

foreach file $files {
   mol load pdb $file
   # Set bead type for end acid after finding COG
    set terminal_res1 [atomselect top "resid 1 and name C1 O1 O2"]
    set com1 [measure center $terminal_res1] 
    set OM2 [atomselect top "resid 1 and name C1"]
    $OM2 moveto $com1
    $OM2 set name TP2
   
   # Set bead type for end alkane after finding COG
    set terminal_res2 [atomselect top "resid 40 and name C2 C3 C4"]
    set com2 [measure center $terminal_res2] 
    set OM1 [atomselect top "resid 40 and name C2"]
    $OM1 moveto $com2
    $OM1 set name SC1
   
   # Alkane groups
    for {set res 1} {$res < 40} {incr res} {
     set alkane [atomselect top "name C2 C3 C4 and resid $res"]
     set com_a [measure center $alkane]
     set AM [atomselect top "name C3 and resid $res"]
     $AM moveto $com_a
     $AM set name SC2
     }
   # Carboxylate groups
    for {set res 2} {$res < 41} {incr res} {
     set ester [atomselect top "name C1 O1 O2 and resid $res"]
     set com_e [measure center $ester]
     set EM [atomselect top "name C1 and resid $res"]
     $EM moveto $com_e
     $EM set name TN4a
     }
   set sel [atomselect top all]
   set towrite [atomselect top "name TP2 SC1 TN4a SC2"]
   $towrite set resname P3HB
   set new_file [ string trim $file ".pdb"]
   $towrite writepdb $new_file-Martini.pdb
}
exit
