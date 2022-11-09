package require psfgen
package require solvate
topology /home/v0e/ff_all/charmff/toppar_c36_jul20/top_all36_carb.rtf

set length 20 
set name native_cell_dp$length 

#build Center chains
set chain C 
 foreach n  { 0 1 2 4 5 6 7 9 10 11 } {
  resetpsf
  set segment $chain$n
  segment  $segment {
  for {set i 1} {$i<=$length} {incr i} {
   residue $i BGLC 
  }
  }
  for {set i $length} {$i>1} {incr i -1} {
   set j [ expr $i -1]
   patch   14bb $segment:$i $segment:$j
  }
  coordpdb csff-elementary-fibril-$segment.pdb $segment
  guesscoord
  regenerate angles dihedrals
  writepdb $segment.pdb
  writepsf $segment.psf
  unset segment
  resetpsf
 }


#build  Origin chains
set chain O 
 foreach n  { 0 1 2 3 5 6 7 8  } { 
  resetpsf
  set segment $chain$n
  segment  $segment {
  for {set i 1} {$i<=$length} {incr i} {
   residue $i BGLC
  }
  }
  for {set i $length} {$i>1} {incr i -1} {
   set j [ expr $i -1]
   patch   14bb $segment:$i $segment:$j
  }
  coordpdb csff-elementary-fibril-$segment.pdb $segment
  guesscoord
  regenerate angles dihedrals
  writepdb $segment.pdb
  writepsf $segment.psf
  unset segment
  resetpsf
 }



#combine
resetpsf
 foreach i  { 1 2 5 6 7 8 11 12 13} {
readpsf  O$i.psf
coordpdb O$i.pdb
}
 foreach i  { 4 5 6 9 10 11 12 15 16 } {
readpsf  C$i.psf
coordpdb C$i.pdb
}
writepdb $name.pdb	 
writepsf $name.psf
 
#solvate
#solvate $name.psf $name.pdb -t 8 -o $name-solvated 	

exit


