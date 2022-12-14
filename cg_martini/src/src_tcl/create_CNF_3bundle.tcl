package require mol
package require topotools
package require Orient
namespace import Orient::orient


## open files for writing
set outfile1 [open "coords.dat" w]

## get all the Martini mapped files for a given acetylation composition

set dir "."
set files [glob -directory $dir initstruct*Martini.pdb]
set files [lsort -dictionary $files]
puts "******sorted order****"
puts ${files}

### Align principal axes along x,y and z for all molecules (Cellulose looks like a cylindrical cross-section)
foreach file $files {
  set sel [atomselect top all]
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 2] {0 0 1}]
  $sel move $A
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 1] {0 1 0}]
  $sel move $A
  set I [draw principalaxes $sel]
}

### Check for long axis of the molecule and align it with the z axis for all files
foreach file $files {
  mol load pdb $file
  set sel0 [atomselect top all]
  set tcl_precision 2
  set min_max_coord [measure minmax $sel0]
  set dist [vecsub [lindex $min_max_coord 1] [lindex $min_max_coord 0]]

  # First get the max distance cooresponding to cellulose fiber axis

  set max_dist [tcl::mathfunc::max {*}$dist]
  set idx [lsearch $dist $max_dist]
  if { $idx eq 0 } {
     set long { 1 0 0}
  } elseif { $idx eq 1 } {
     set long { 0 1 0}
  } elseif { $idx eq 2 } {
     set long { 0 0 1}
  }
    else {
     puts "index for long axis wrong"
     break;
  }

  # Align longest axis along z axis
  set sel [atomselect top all]
  set M [transvecinv $long]
  $sel move $M
  set M [transaxis y -90]
  $sel move $M

}

## Set all the files to center fiber at the origin
set file_counter 1
foreach file $files {
  mol load pdb $file
  set A [atomselect top all]
  set minus_com [vecsub {0.0 0.0 0.0} [measure center $A]]
  $A moveby $minus_com
  set new_file "initstruct-s$file_counter-Martini-centered.pdb"
  $A writepdb $new_file
  incr file_counter 1
}

set dir "."
set files [glob -directory $dir initstruct*centered.pdb]
set files [lsort -dictionary $files]
set file_counter 0

foreach file $files {
  incr file_counter 1
  mol load pdb $file
  set sel0 [atomselect top all]
  set tcl_precision 2
  set min_max_coord [measure minmax $sel0]
  set dist [vecsub [lindex $min_max_coord 1] [lindex $min_max_coord 0]]
  set max_dist [tcl::mathfunc::max {*}$dist]
  set idx [lsearch $dist $max_dist]
  # Delete longest dimension in the distance list corresponding to fiber axis
  set dist [lreplace $dist $idx $idx] 
  # Find the maximum among short axis to place the CNFs in a regular geometry
  set max_dist [tcl::mathfunc::max {*}$dist]
  puts $outfile1 "$min_max_coord\ $dist\ $max_dist"
  set tcl_precision 2
  set fiber_spacing 10.0
  # set moveby coordinates to move CNFs relative to the one in the center
  set d [expr $max_dist + $fiber_spacing]
  set new_file "initstruct-s$file_counter-Martini-moved.pdb"
  # write the moved coordinates to a different pdb
  if {$file_counter eq 1} {
     ### write dummy file for format sake, not actually moved since this is the central fiber
     set sel [atomselect top all]
     $sel writepdb $new_file
     puts "success"
     continue
  } elseif {$file_counter eq 2} {
    set sel [atomselect top all]
    $sel moveby [list [expr {$d * 1}] 0 0]
  } elseif {$file_counter eq 3} {
    set sel [atomselect top all]
    $sel moveby [list [expr {$d * 0.5 * -1}] [expr {$d * 0.5}] 0]
  }
  $sel writepdb $new_file
}

set dir "."
set files [glob -directory $dir initstruct*moved.pdb]
set files [lsort -dictionary $files]
set file_counter 0

foreach file $files {
  mol load pdb $file
  set sel_$file_counter [atomselect top all]
  incr file_counter 1
  }
set final_sel [::TopoTools::selections2mol "$sel_0 $sel_1 $sel_2"]

animate write pdb acetylated_CNF_bundle.pdb $final_sel

# final check to make sure it is centered at the origin
mol load pdb acetylated_CNF_bundle.pdb
set A [atomselect top all]
set minus_com [vecsub {0.0 0.0 0.0} [measure center $A]]
$A moveby $minus_com
$A writepdb acetylated_CNF_bundle.pdb
$A writegro acetylated_CNF_bundle.gro

exit
