package require mol

set dir [lindex $argv 0]
puts "$dir"
set files [glob -directory "$dir" *.pdb]
set files [lsort -dictionary $files
]
foreach file $files {
    mol load pdb $file
   # Set bead type for end acid after finding COM
   
    set terminal_res1 [atomselect top "resid 1 and name C1 O1 O2"]
    set com1 [measure center $terminal_res1] 
    set OM2 [atomselect top "resid 1 and name C1"]
    $OM2 moveto $com1
    $OM2 set name TP2

   # Set bead type for end alcohol after finding COM
   
    set terminal_res2 [atomselect top "resid 20 and name C10 OZ"]
    set com2 [measure center $terminal_res2] 
    set OM1 [atomselect top "resid 20 and name C10"]
    $OM1 moveto $com2
    $OM1 set name TP1
   
   # aromatic atoms
    for {set res 1} {$res < 21} {incr res 1} {
     # first set
     set ring_1 [atomselect top "name C2 C5 and resid $res"]
     set com_r1 [measure center $ring_1]
     set RM1 [atomselect top "name C2 and resid $res"]
     $RM1 moveto $com_r1
     $RM1 set name TC5
     # second set
     set ring_2 [atomselect top "name C3 C4 and resid $res"]
     set com_r2 [measure center $ring_2]
     set RM2 [atomselect top "name C3 and resid $res"]
     $RM2 moveto $com_r2
     $RM2 set name TC5
    # third set
     set ring_3 [atomselect top "name C6 C7 and resid $res"]
     set com_r3 [measure center $ring_3]
     set RM3 [atomselect top "name C6 and resid $res"]
     $RM3 moveto $com_r3
     $RM3 set name TC5
     }
   # Ester groups of one residue
    for {set res 1} {$res < 20} {incr res 1} {
     set ester_1 [atomselect top "name C8 C9 O3 O4 and resid $res"]
     set com_e1 [measure center $ester_1]
     set EM1 [atomselect top "name C8 and resid $res"]
     $EM1 moveto $com_e1
     $EM1 set name N4a
     }
   # Ester groups of adjoining residues
    for {set i 1; set j 2} {$i < 20} {incr i; incr j} {
     set ester_2 [atomselect top "(name C10 and resid $i) or (name O1 C1 O2 and resid $j)"]
     set com_e2 [measure center $ester_2]
     set EM2 [atomselect top "name C10 and resid $i"]
     $EM2 moveto $com_e2
     $EM2 set name N4a
     }
   # Ester groups of last residue
    set ester_3 [atomselect top "name C8 C9 O3 O4 and resid 20"]
    set com_e3 [measure center $ester_3]
    set EM3 [atomselect top "name C8 and resid 20"]
    $EM3 moveto $com_e3
    $EM3 set name N4a
   
   set towrite [atomselect top "name TP2 TP1 N4a TC5"]
   $towrite set resname PETG 
   set new_file [ string trim $file ".pdb"]
   $towrite writepdb $new_file-Martini.pdb
   $towrite writegro $new_file-Martini.gro
}
exit
