# Generate input file for packmol to pack CNF bundle and polymer

## Set the paths to the files
set dir1 [lindex $argv 0]
set dir2 [lindex $argv 1]
set nchains [lindex $argv 2]

## get number of polymer files and divide number of chains to have packmol inputs 
set polymer_files [glob -directory "$dir2" *.pdb]
set polymer_files [lsort -dictionary $polymer_files]
set n_poly_files [llength $polymer_files]
set n_per_file [expr {$nchains/$n_poly_files}]
set n_per_file [expr {int($n_per_file)}]
set n_last_file [expr {$nchains - $n_per_file * [expr {$n_poly_files - 1}]}]

puts "**********\ $n_poly_files********************"

set cnf_file [glob -directory "$dir1" acetylated_CNF_bundle.pdb]
set cnf_file [lsort -dictionary $cnf_file]

# Find maximum and minimum coordinates of the bundle to set box size

mol load pdb $cnf_file
set sel [atomselect top all]
set min_max_coord [measure minmax $sel]
set dist [vecsub [lindex $min_max_coord 1] [lindex $min_max_coord 0]]
#puts $outfile1 "bundle:\ $min_max_coord\ $dist"

set xlo [lindex $min_max_coord 0 0]
set ylo [lindex $min_max_coord 0 1]
set zlo [lindex $min_max_coord 0 2]
set xhi [lindex $min_max_coord 1 0]
set yhi [lindex $min_max_coord 1 1]
set zhi [lindex $min_max_coord 1 2]

set xdist [lindex $dist 0]
set ydist [lindex $dist 1]
set zdist [ lindex $dist 2]

set newxdist [ expr { $xdist + 0.1*$xdist } ]
set newydist [ expr { $ydist + 0.1*$ydist } ]
set newzdist [ expr { $zdist + 0.1*$zdist } ]

set xlo [ expr {$xlo - $newxdist/2} ]
set ylo [ expr {$ylo - $newydist/2} ]
set zlo [ expr {$zlo - $newzdist/2} ]

set xhi [ expr {$xhi + $newxdist/2} ]
set yhi [ expr {$yhi + $newydist/2} ]
set zhi [ expr {$zhi + $newzdist/2} ]


# write packmol file
set pack_file [open "packmol_CNF_polymer.inp" w]
puts $pack_file "tolerance\ 2"
puts $pack_file "filetype\ pdb"
puts $pack_file "output\ acetylated_CNF_polymer.pdb"
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
