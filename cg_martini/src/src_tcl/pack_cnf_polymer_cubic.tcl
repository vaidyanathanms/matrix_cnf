# Generate input file for packmol to pack CNF bundle and polymer

## Set the paths to the files
set dir1 [lindex $argv 0]
set dir2 [lindex $argv 1]
set nchains [lindex $argv 2]
set ncnf_per_bundle [lindex $argv 3]

## get number of polymer files and divide number of chains to have packmol inputs 
set polymer_files [glob -directory "$dir2" *.pdb]
set polymer_files [lsort -dictionary $polymer_files]
set n_poly_files [llength $polymer_files]
set n_per_file [expr {$nchains/$n_poly_files}]
set n_per_file [expr {int($n_per_file)}]
set n_last_file [expr {$nchains - $n_per_file * [expr {$n_poly_files - 1}]}]

puts "**********\ $n_poly_files********************"

set cnf_file [glob -directory "$dir1" *.pdb]
set cnf_file [lsort -dictionary $cnf_file]

# Find maximum and minimum coordinates of the bundle to set box size

mol load pdb $cnf_file
set sel [atomselect top all]
set min_max_coord [measure minmax $sel]
set dist [vecsub [lindex $min_max_coord 1] [lindex $min_max_coord 0]]
set max_dist [tcl::mathfunc::max {*}$dist]
set buffer_dist 20.0
set new_dist [expr {$max_dist/2 + $buffer_dist/2} ]
#puts $outfile1 "bundle:\ $min_max_coord\ $dist"

set xlo [ expr { (-1 * $new_dist) - ($buffer_dist * ($ncnf_per_bundle - 1) ) } ]
set ylo [ expr { (-1 * $new_dist) - ($buffer_dist * ($ncnf_per_bundle - 1) ) }]
set zlo [ expr { (-1 * $new_dist) - ($buffer_dist * ($ncnf_per_bundle - 1) ) }]

set xhi [ expr { $new_dist + ($buffer_dist * ($ncnf_per_bundle - 1)) } ]
set yhi [ expr { $new_dist + ($buffer_dist * ($ncnf_per_bundle - 1)) } ]
set zhi [ expr { $new_dist + ($buffer_dist * ($ncnf_per_bundle - 1)) } ]

# write packmol file
set pack_file [open "packmol_CNF_polymer.inp" w]
puts $pack_file "tolerance\ 2"
puts $pack_file "filetype\ pdb"
puts $pack_file "output\ CNF_polymer.pdb"
puts $pack_file "seed\ -1"
puts $pack_file "\n"
puts $pack_file "#Begin\ structure\ generation"
puts $pack_file "structure\ $cnf_file"
puts $pack_file "\  number\   1"
puts $pack_file "\  center"
puts $pack_file "\ fixed\ 0\  0\  0\  0\  0\  0"
puts $pack_file "end\ structure"
puts $pack_file "\n"

set ctr 0
set last_file_ctr [expr {$n_poly_files - 1}] 
foreach file $polymer_files {
  if {$ctr eq $last_file_ctr} {
      puts $pack_file "structure\  $file"
      puts $pack_file "\  number\   $n_last_file"
      puts $pack_file "\ inside\ box\ $xlo\  $ylo\  $zlo\  $xhi\  $yhi\  $zhi"
      puts $pack_file "end\ structure"
      puts $pack_file "\n"
  } else {
  puts $pack_file "structure\  $file"
  puts $pack_file "\  number\   $n_per_file"
  puts $pack_file "\ inside\ box\ $xlo\  $ylo\  $zlo\  $xhi\  $yhi\  $zhi"
  puts $pack_file "end\ structure"
  puts $pack_file "\n"
  }
  incr ctr 1
}
exit
#
