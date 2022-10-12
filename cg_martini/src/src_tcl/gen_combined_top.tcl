### To generate combined topology file for cnf polymer systems

## Set the paths to the files
set dir1 [lindex $argv 0]
set dir2 [lindex $argv 1]
set restype [lindex $argv 2] 
set nchains [lindex $argv 3]
set ncnf_per_bundle [lindex $argv 4]

set cnf_topfiles [glob -directory "$dir1" initstruct-*-acetylated.top]
set cnf_topfiles [lsort -dictionary $cnf_topfiles]
set polymer_topfiles [glob -directory "$dir2" *.itp] 
set polymer_topfiles [lsort -dictionary $polymer_topfiles]

set outfile1 [open "alltop.top" w]

puts $outfile1 ";;\ -------------------------------------"
puts $outfile1 ";;\ Combined topology file"
puts $outfile1 ";;\ Generated\ using\ create_packmol_file.py"
puts $outfile1 ";;\ Use\ with\ GROMACS\ grompp"
puts $outfile1 ";;\ -------------------------------------"
puts $outfile1 ";;\ Include\ interaction\ parameter\ files"

puts $outfile1 "#include\ \"/lustre/or-scratch/cades-bsd/world-shared/cnf_martini_cg/ff_files/cell_martini/cell_martini3.itp\""
foreach fname ${cnf_topfiles} {
 puts $outfile1 "#include\ \"${fname}\""
}

foreach fname ${polymer_topfiles} {
 puts $outfile1 "#include\ \"${fname}\""
} 

puts $outfile1 "\n"
puts $outfile1 "\[\ System\ \]"
puts $outfile1 ";\ Name"
puts $outfile1 "\ Combined"
puts $outfile1 "\n"
puts $outfile1 "\[\ molecules\ \]"
puts $outfile1 ";\ Compound\      #mols"
for {set i 1} {$i <= $ncnf_per_bundle} {incr i 1} {
    puts $outfile1 "CELLULOSE-acetylated-$i\    1"
}
puts $outfile1 ";matrix"
puts $outfile1 "$restype\    $nchains"

exit
