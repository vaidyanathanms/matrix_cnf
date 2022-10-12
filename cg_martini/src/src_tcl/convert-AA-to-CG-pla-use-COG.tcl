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
   
   # Set bead type for end alcohol after finding COG
   
    set terminal_res2 [atomselect top "resid 40 and name C2 C3"]
    set com2 [measure center $terminal_res2] 
    set AM1 [atomselect top "resid 40 and name C2"]
    $AM1 moveto $com2
    $AM1 set name TC1
   
   # Ester groups of adjoining residues
    for {set i 1; set j 2} {$i < 40} {incr i; incr j} {
     set ester [atomselect top "(name C2 C3 and resid $i) or (name O1 C1 O2 and resid $j)"]
     set com_e [measure center $ester]
     set EM [atomselect top "name C3 and resid $i"]
     $EM moveto $com_e
     $EM set name N4a
     }
   set sel [atomselect top all]
   #$sel writepdb test.pdb
   set towrite [atomselect top "name TP2 TC1 N4a"]
   $towrite set resname PLA
   set new_file [ string trim $file ".pdb"]
   $towrite writepdb $new_file-Martini.pdb
   $towrite writegro $new_file-Martini.gro
}
exit
