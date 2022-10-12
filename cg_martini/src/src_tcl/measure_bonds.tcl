
mol load initconf.gro
set outfile1 [open "bond_dist.dat" w]
set sel [atomselect top all]
set natoms [eval $sel num]
for {set i 0; set j 1} {$atid < [expr $natoms - 1]} {incr i 1; incr j 1} {
    set bnd_l [measure bond $i $j]
    puts outfile1 $bnd_l
}
