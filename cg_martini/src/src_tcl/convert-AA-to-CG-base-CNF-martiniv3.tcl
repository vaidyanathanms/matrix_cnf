## map atomistic structure to Martini v3.0 beads for umodified Cellulose fibers##
### command: vmd -dispdev text -e convert-AA-to-CG-base-CNF-martiniv3.tcl <file_type> <coord_file_name> <nchains> <DP> <cnfs_per_bundle>
### file_type - pdb or gro
### coord_file_name - input PDB/GRO file
### nchains - number of cellulose chains per CNF (18, 36, etc)
### DP - number of glycan monomers in a cellulose strand
### cnfs_per_bundle - number of CNFs in the bundle you would like to generate
##------------------------------------------------------------------------ 
package require mol
######### Read the input coordinate file in pdb/gro format
set file_type [lindex $argv 0]
set coord_file [lindex $argv 1]

if { $file_type eq "pdb" } {
   mol load pdb $coord_file 
} elseif { $file_type eq "gro" } {
   mol load gro $coord_file
} else {
   puts "incorrect file type"
   exit
}
######## Set the system information
set nchains [lindex $argv 2]
set DP [lindex $argv 3]
# for any given monomer without acetylation, the maximum atoms per resid cannot exceed 23, when it does check for acetylated atom types 
set max_atoms_per_res 23
set cnfs_per_bundle [lindex $argv 4]
set nfragments [expr $nchains * $cnfs_per_bundle]
set nsegments [expr $nchains * $cnfs_per_bundle]

# Bead types 
# B1 (hydroxymethyl group)
# B2 (ring oxygen and anomeric carbon) if it is not a reducing end, include hydroxyl group for the a reducing end 
# B3 (vicinal diol)
# B4 Hydroxyethyl group (include hydroxyl hydrogen for the non-reducing end)
###### Note: In this code, input structure is assumed to have reducing end at resid 1 and non-reducing end at resid 20

for {set frag 0} {$frag < $nfragments} {incr frag 1} {
   for {set res 1} {$res <= $DP} {incr res 1} {
      set b1 [atomselect top "(name O6 HO6 C6 H61 H62) and fragment $frag and resid $res"]
      set com1 [measure center $b1 weight mass]
      set c6 [atomselect top "name C6 and fragment $frag and resid $res"]
      $c6 moveto $com1
      $c6 set name B1
      if {$res eq 1} {
         set b2 [atomselect top "(name O5 C1 H1 O1 HO1) and fragment $frag and resid $res"]
         set com2 [measure center $b2 weight mass]
         set o5 [atomselect top "name O5 and fragment $frag and resid $res"]
         $o5 moveto $com2
         $o5 set name B2
       } else {
         set b2 [atomselect top "(name O5 C1 H1) and fragment $frag and resid $res"]
         set com2 [measure center $b2 weight mass]
         set o5 [atomselect top "name O5 and fragment $frag and resid $res"]
         $o5 moveto $com2
         $o5 set name B2
       }
      set b3 [atomselect top "(name C2 H2 C3 H3 O2 HO2 O3 HO3) and fragment $frag and resid $res"]
      set com3 [measure center $b3 weight mass]
      set o2 [atomselect top "name O2 and fragment $frag and resid $res"]
      $o2 moveto $com3 
      $o2 set name B3
      if {$res eq $DP} {
         set b4 [atomselect top "(name C5 H5 C4 H4 O4 HO4) and fragment $frag and resid $res"]
         set com4 [measure center $b4 weight mass]
         set o4 [atomselect top "name O4 and fragment $frag and resid $res"]
         $o4 moveto $com4
         $o4 set name B4
       } else {
         set b4 [atomselect top "(name C5 H5 C4 H4 O4) and fragment $frag and resid $res"]
         set com4 [measure center $b4 weight mass]
         set o4 [atomselect top "name O4 and fragment $frag and resid $res"]
         $o4 moveto $com4
         $o4 set name B4
       }
   }
}
set T1 [atomselect top "name B1"]
$T1 set name T1
set T2 [atomselect top "name B2"]
$T2 set name T2
set R3 [atomselect top "name B3"]
$R3 set name R3
set S4 [atomselect top "name B4"]
$S4 set name S4

set sel_all [atomselect top "name T1 T2 R3 S4"]
$sel_all set resname CELL

$sel_all writepdb base-CNF-bundle.pdb

#### Change the order of atoms consistent with the top file (GLCB14.top) obtained from carbo2martini3.py code 

mol load pdb base-CNF-bundle.pdb
for {set seg 1} {$seg <= $nsegments} {incr seg 1} {
   for {set res 1} {$res <= $DP} {incr res 1} {
      set S4 [atomselect top "name S4 and segname C$seg and resid $res"]
      set coord_S4 [$S4 get {x y z}]
      set T1 [atomselect top "name T1 and segname C$seg and resid $res"]
      set coord_T1 [$T1 get {x y z}]
      $S4 moveto [expr $coord_T1]
      $T1 moveto [expr $coord_S4]
      $T1 set name S4
      $S4 set name T1
      set R3 [atomselect top "name R3 and segname C$seg and resid $res"]
      set coord_R3 [$R3 get {x y z}]
      set T1 [atomselect top "name T1 and segname C$seg and resid $res"]
      set coord_T1 [$T1 get {x y z}]
      $R3 moveto [expr $coord_T1]
      $T1 moveto [expr $coord_R3]
      $T1 set name R3
      $R3 set name T1
      set T2 [atomselect top "name T2 and segname C$seg and resid $res"]
      set coord_T2 [$T2 get {x y z}]
      set T1 [atomselect top "name T1 and segname C$seg and resid $res"]
      set coord_T1 [$T1 get {x y z}]
      $T2 moveto [expr $coord_T1]
      $T1 moveto [expr $coord_T2]
      $T1 set name T2
      $T2 set name T1

   }
}

set sel [atomselect top "resid 20 and name T2"]
$sel set name S2
set sel [atomselect top all]
$sel writepdb base-CNF-bundle.pdb
$sel writegro base-CNF-bundle.gro

exit
