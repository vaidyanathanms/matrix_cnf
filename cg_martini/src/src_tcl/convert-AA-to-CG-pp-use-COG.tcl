
package require mol

set dir [lindex $argv 0]
puts "$dir" 
set files [glob -directory "$dir" *.pdb]
set files [lsort -dictionary $files]

foreach file $files {
    mol load pdb $file
    # Set bead type for terminal propane
    set terminal_res1 [atomselect top "resid 1 and name C1 C2 C3"]
    set com1 [measure center $terminal_res1] 
    set C_end1 [atomselect top "resid 1 and name C1"]
    $C_end1 moveto $com1
    $C_end1 set name SC1

    # Set bead type for terminal isopropanol
    set terminal_res2 [atomselect top "resid 40 and name C1 C2 C3"]
    set com2 [measure center $terminal_res2] 
    set C_end2 [atomselect top "resid 40 and name C1"]
    $C_end2 moveto $com2
    $C_end2 set name SC1

    # Set bead type for monomers (middle)
    for {set res 2} {$res < 40} {incr res 1} {
    set middle_res [atomselect top "resid $res and name C1 C2 C3"]
    set com_3 [measure center $middle_res]
    set C_middle [atomselect top "resid $res and name C1"]
    $C_middle moveto $com_3
    $C_middle set name SC1
    }
    set sel_all [atomselect top "name SC1"]
    $sel_all set resname PP 
    set new_file [ string trim $file ".pdb"]
    $sel_all writepdb $new_file-Martini.pdb
    $sel_all writegro $new_file-Martini.gro
}
exit
