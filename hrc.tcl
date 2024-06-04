# used in VMD https://www.ks.uiuc.edu/Research/vmd/
# 
# calculate bond length and ring area distributions
# file_1 file_2: file names of the trajectory, for most files, file_1 = file_2, 
#               for psf and dcd paird files, file_1 = psf file, file_2 = dcd file               
# mode: bond, calculate bond length only; ring, calculate ring area only.
# firstFrame lastFrame step: user-specified trajectory part
# cutoff: the cutoff used to determine a bond
# selection: the atoms that will be involved in the calculation
# example: source hrc.tcl; hrc test.pdb test.pdb ring 0 9 1 1.6 {name C}

# NOTE: for the bond mode, it calculates all bonds information in the selection, not only these appear in rings.

proc hrc {file_1 file_2 mode firstFrame lastFrame step cutoff selection} {

if {[string equal $file_1 $file_2]} {
	mol new $file_1 first $firstFrame last $lastFrame step $step waitfor all
} else {
	mol new $file_1
	mol addfile $file_2 first $firstFrame last $lastFrame step $step waitfor all
}

set logfile [open log.dat w]
puts $logfile "********************** file: $file_1 $file_2 **********************"
puts $logfile "Command: hrc $file_1 $file_2 $mode $firstFrame $lastFrame $step $cutoff {$selection}"

set frameNum [molinfo top get numframes]

puts "   "
puts "frameNum = $frameNum"
puts "frameID  number  avgRingAreaSingle  avgRingCircularitySingle"

puts $logfile "frameNum= $frameNum"
puts $logfile "   "

set bondFileStyle 1
set ringFileStyle 1

# frame loop of the file
for {set n 0} {$n < $frameNum} {incr n} {
	set frameID [expr {$n*$step+$firstFrame}]
	set atomsel [atomselect top "$selection" frame $n]
	set atomList [$atomsel get index]
	set atomNum [llength $atomList]
	set bondLengthRecordSingle {}
	set ringAreaRecordSingle {}
	
	# update bond connection according to user-specified cutoff-----------------------------
	topo clearbonds
	foreach element $atomList {
		set bonded_atoms_sel [atomselect top "$selection and within $cutoff of index $element" frame $n]
		set bonded_atoms_list [$bonded_atoms_sel list]
		foreach bonded_atom $bonded_atoms_list {
			if {$bonded_atom != $element} {
				topo addbond $element $bonded_atom
			}
		}	
	}
	#---------------------------------------------------------------------------------------
		
	set bondList [topo getbondlist -sel $atomsel]
	set bondNumSingle [llength $bondList]
	puts $logfile "   "
	puts $logfile "---------- frameID= $frameID, atomNum= $atomNum, bondNumSingle= $bondNumSingle ----------"
	puts $logfile "bondList = "
	puts $logfile "$bondList"
	
	# bond mode --------------------------------------------------------------------------------------- start
	if {[string equal -nocase $mode bond]} {
		# set file writing style
		if {$bondFileStyle == 1} {
			set bondLengthAvgFile [open bondLengthAvgFile.dat w]
			puts $bondLengthAvgFile "#frameID  bondNumSingle  avgBondLengthSingle"
			set bondLengthAllFile [open bondLengthAllFile.dat w]
			set bondFileStyle 0
			close $bondLengthAvgFile	
			close $bondLengthAllFile
		}
		set bondLengthAvgFile [open bondLengthAvgFile.dat a]
		set bondLengthAllFile [open bondLengthAllFile.dat a]
		
		for {set i 0} {$i < $bondNumSingle} {incr i} {
			set bondLength [measure bond [lindex $bondList $i] first $n last $n]
			lappend bondLengthRecordSingle $bondLength
			puts $bondLengthAllFile [format "%-10.4f" $bondLength] ; # write bond length to file
		}
		
		set bondLengthSumSingle 0
		foreach element $bondLengthRecordSingle {
			set bondLengthSumSingle [expr {$bondLengthSumSingle + $element}]
		}
		
		set avgBondLengthSingle [expr {$bondLengthSumSingle/$bondNumSingle}]
		puts $bondLengthAvgFile [format "%-10d %-10d %-10.4f" $frameID  $bondNumSingle  $avgBondLengthSingle]
		puts [format "%-10d %-10d %-10.4f" $frameID $bondNumSingle $avgBondLengthSingle]
		close $bondLengthAvgFile	
		close $bondLengthAllFile
	# bond mode --------------------------------------------------------------------------------------- end
	
	# ring mode --------------------------------------------------------------------------------------- start	
	} elseif {[string equal -nocase $mode ring]} {
		# set file writing style
		if {$ringFileStyle == 1} {
			set ringAvgFile [open ringAvgFile.dat w]
			puts $ringAvgFile "#frameID  ringNumSingle  avgRingAreaSingle  avgRingCircularitySingle"
			set ringAllFile [open ringAllFile.dat w]
			set ringFileStyle 0
			close $ringAvgFile	
			close $ringAllFile
		}
		set ringAvgFile [open ringAvgFile.dat a]
		set ringAllFile [open ringAllFile.dat a]

		set ringNumSingle 0
		set ringListSingle {}
		set ringListSortedSingle {}
		
		# -------------------------------------- start to search the rings
		foreach element $atomList {
			# -------------------------------------- assume the 1st atom of the ring
			set atom_1_sel [atomselect top "index $element" frame $n]
			set x1 [$atom_1_sel get x]
			set y1 [$atom_1_sel get y]
			set z1 [$atom_1_sel get z]
			
			# -------------------------------------- search the 2nd and 3rd atoms of the ring
			set bonded_atoms_list [$atom_1_sel getbonds]
			set bonded_atoms_list [lindex $bonded_atoms_list 0]
			
			set count [llength $bonded_atoms_list]
			if {$count <= 1} {
				continue
				# If count==1, it means this situation cannot form a ring;
				# only if count==2 or ==3, it has chance to form a ring
			} elseif {$count == 2 || $count == 3} {
				set aa [lindex $bonded_atoms_list 0]
				set bb [lindex $bonded_atoms_list 1]
				set cc [lindex $bonded_atoms_list 2]
				set situationList [list [list $aa $bb] [list $bb $cc] [list $aa $cc]]
				# In this foreach loop, the first atom may bond with 2 or 3 other atoms, so it will have 3 situations;
				# if there is no cc (count==2), one situation will have two index values, and another
				# two situations will have an index value and an empty value;
				# so, for no cc condition, just find and remove the two situations by if continue commands, 
				# and only consider the left situation which has two index values.
				# If count==3, all the three situations will be considered by command foreach.
				# 这确实是个比较烧脑的逻辑 !>_<
				foreach situation $situationList {
					set ringDone 0
					if {[llength [lindex $situation 0]] == 0 || [llength [lindex $situation 1]] == 0} {
						continue
					}
					set atom_2_id [lindex $situation 0]
					set atom_3_id [lindex $situation 1]
					set atom_2_sel [atomselect top "index $atom_2_id" frame $n]
					set atom_3_sel [atomselect top "index $atom_3_id" frame $n]
					set x2 [$atom_2_sel get x]
					set y2 [$atom_2_sel get y]
					set z2 [$atom_2_sel get z]
					set x3 [$atom_3_sel get x]
					set y3 [$atom_3_sel get y]
					set z3 [$atom_3_sel get z]
					
					# -------------------------------------- search the 4th atom of the ring
					set bonded_atoms_list [$atom_2_sel getbonds]
					set bonded_atoms_list [lindex $bonded_atoms_list 0]
					# "self" is the previous atom, in this context, the first atom which is bonded with atom_2 and atom_3,
					# because we want to find who are bonded to atom_2, atom_1 must be one, so name it with "self".
					set self [lsearch $bonded_atoms_list $element]
					# In "bonded_atoms_list_new", we delete self, so the rest are all new atoms.
					set bonded_atoms_list_new [lreplace $bonded_atoms_list $self $self]			
					set count [llength $bonded_atoms_list_new ]
					if {$count == 0} {
						continue
					} elseif {$count == 1} {
						set atom_4_id [lindex $bonded_atoms_list_new 0]
						set atom_4_sel [atomselect top "index $atom_4_id" frame $n]
						set x4 [$atom_4_sel get x]
						set y4 [$atom_4_sel get y]
						set z4 [$atom_4_sel get z]
						# Using the angle to determine is atom_4 is within the ring
						set angle_ref [angle $x1 $y1 $z1 $x3 $y3 $z3 $x2 $y2 $z2]
						set angle_4 [angle $x1 $y1 $z1 $x3 $y3 $z3 $x4 $y4 $z4]
						if {$angle_4 > $angle_ref} {
							continue
						}	
					} elseif {$count == 2} {
						set atom_40_id [lindex $bonded_atoms_list_new 0]
						set atom_40_sel [atomselect top "index $atom_40_id" frame $n]
						set x40 [$atom_40_sel get x]
						set y40 [$atom_40_sel get y]
						set z40 [$atom_40_sel get z]
						set atom_41_id [lindex $bonded_atoms_list_new 1]
						set atom_41_sel [atomselect top "index $atom_41_id" frame $n]
						set x41 [$atom_41_sel get x]
						set y41 [$atom_41_sel get y]
						set z41 [$atom_41_sel get z]
						
						set angle_40 [angle $x1 $y1 $z1 $x3 $y3 $z3 $x40 $y40 $z40]
						set angle_41 [angle $x1 $y1 $z1 $x3 $y3 $z3 $x41 $y41 $z41]
						# for count==2 condition, there must be one atom the ring atom and another one not
						if {$angle_40 < $angle_41} {
							set atom_4_id $atom_40_id
						} else {
							set atom_4_id $atom_41_id	
						}
						set atom_4_sel [atomselect top "index $atom_4_id" frame $n] 
						set x4 [$atom_4_sel get x]
						set y4 [$atom_4_sel get y]
						set z4 [$atom_4_sel get z]
					}
					
					# --------------------------------------find the 5th atom of the ring
					# Same precedure with atom_4.
					set bonded_atoms_list [$atom_3_sel getbonds]
					set bonded_atoms_list [lindex $bonded_atoms_list 0]		
					set self [lsearch $bonded_atoms_list $element]
					set bonded_atoms_list_new [lreplace $bonded_atoms_list $self $self]			
					set count [llength $bonded_atoms_list_new ]
					if {$count == 0} {
						continue
					} elseif {$count == 1} {
						set atom_5_id [lindex $bonded_atoms_list_new 0]
						set atom_5_sel [atomselect top "index $atom_5_id" frame $n]
						set x5 [$atom_5_sel get x]
						set y5 [$atom_5_sel get y]
						set z5 [$atom_5_sel get z]
						set angle_ref [angle $x1 $y1 $z1 $x3 $y3 $z3 $x2 $y2 $z2]
						set angle_5 [angle $x1 $y1 $z1 $x2 $y2 $z2 $x5 $y5 $z5]
						if {$angle_5 > $angle_ref} {
							continue
						}	
					} elseif {$count == 2} {
						set atom_50_id [lindex $bonded_atoms_list_new 0]
						set atom_50_sel [atomselect top "index $atom_50_id" frame $n]
						set x50 [$atom_50_sel get x]
						set y50 [$atom_50_sel get y]
						set z50 [$atom_50_sel get z]
						set atom_51_id [lindex $bonded_atoms_list_new 1]
						set atom_51_sel [atomselect top "index $atom_51_id" frame $n]
						set x51 [$atom_51_sel get x]
						set y51 [$atom_51_sel get y]
						set z51 [$atom_51_sel get z]
						
						set angle_50 [angle $x1 $y1 $z1 $x2 $y2 $z2 $x50 $y50 $z50]
						set angle_51 [angle $x1 $y1 $z1 $x2 $y2 $z2 $x51 $y51 $z51]
						if {$angle_50 < $angle_51} {
							set atom_5_id $atom_50_id
						} else {
							set atom_5_id $atom_51_id
						}
						set atom_5_sel [atomselect top "index $atom_5_id" frame $n] 
						set x5 [$atom_5_sel get x]
						set y5 [$atom_5_sel get y]
						set z5 [$atom_5_sel get z]
					}
				
					# -------------------------------------- search the 6th atom of a ring
					set bonded_atoms_list_4 [$atom_4_sel getbonds]
					set bonded_atoms_list_4 [lindex $bonded_atoms_list_4 0]
					set self [lsearch $bonded_atoms_list_4 $atom_2_id]
					set bonded_atoms_list_4_new [lreplace $bonded_atoms_list_4 $self $self]
					set count_4 [llength $bonded_atoms_list_4_new]
					
					set bonded_atoms_list_5 [$atom_5_sel getbonds]
					set bonded_atoms_list_5 [lindex $bonded_atoms_list_5 0]
					set self [lsearch $bonded_atoms_list_5 $atom_3_id]
					set bonded_atoms_list_5_new [lreplace $bonded_atoms_list_5 $self $self]
					set count_5 [llength $bonded_atoms_list_5_new]
					
					# 
					if {$count_4>0 && $count_5>0} {
						set FLAG 0
						foreach index_4 $bonded_atoms_list_4_new {
							if {$FLAG == 1} {
								break
							} else {
								foreach index_5 $bonded_atoms_list_5_new {
									 if {$index_4 == $index_5} {
										set FLAG 1
										set ringDone 1
										# When FLAG==1, index_4 or index_5 is atom_6.
										set atom_6_id $index_4
										set atom_6_sel [atomselect top "index $atom_6_id" frame $n] 
										set x6 [$atom_6_sel get x]
										set y6 [$atom_6_sel get y]
										set z6 [$atom_6_sel get z]
										break
									} else {
										continue
									}
								}; # foreach index_5 $bonded_atoms_list_5_new					
							}
						}; # foreach index_4 $bonded_atoms_list_4_new
					} else {
						continue; # foreach situation $situationList; continue to next situation loop	
					}
					# 
					
					if {$ringDone == 1} {
						# NOTE: in the same frame, when loop all atoms and search rings, there would be many overlaps
						# ringListSingle: store all found ring_atom_index in the original atom order, even overlap
						# note: the indexes are in ring-direction sequence
						lappend ringListSingle [list $element $atom_2_id $atom_4_id $atom_6_id $atom_5_id $atom_3_id]
						# ringListSingle: store all found ring_atom_index in the sorted atom order, even overlap
						lappend ringListSortedSingle [lsort -integer [list $element $atom_2_id $atom_4_id $atom_6_id $atom_5_id $atom_3_id]]
					}

				} ; # foreach situation $situationList				
			} ; #  elseif {$count == 2 || $count == 3}
		} ; # atom loop
		
# -----------------------------------------------------------------------------------------------------------------------20240320		
		set uniqueRingListSingle {}		; # store the unique ring_atom_index in a single frame, taken from the sorted list
		set uniqueRingAreaSingle {}	; # store the unique ring areas of a frame
		set uniqueRingCircularitySingle {}	; # store the unique ring circularities of a frame
		foreach item $ringListSortedSingle {
			if {$item ni $uniqueRingListSingle} {
				lappend uniqueRingListSingle $item
				set uniqueRingIndex [lsearch $ringListSortedSingle $item]
				# NOTE: the indexes of each atom below is DIFFERENT from the above,
				#       they are atom indexes from the original element atom with a counter-clockwise direction
				set atom_1_id [lindex $ringListSingle $uniqueRingIndex 0]
				set atom_1_sel [atomselect top "index $atom_1_id" frame $n]
				set x1 [$atom_1_sel get x]
				set y1 [$atom_1_sel get y]
				set z1 [$atom_1_sel get z]
				set atom_2_id [lindex $ringListSingle $uniqueRingIndex 1]
				set atom_2_sel [atomselect top "index $atom_2_id" frame $n]
				set x2 [$atom_2_sel get x]
				set y2 [$atom_2_sel get y]
				set z2 [$atom_2_sel get z]
				set atom_3_id [lindex $ringListSingle $uniqueRingIndex 2]
				set atom_3_sel [atomselect top "index $atom_3_id" frame $n]
				set x3 [$atom_3_sel get x]
				set y3 [$atom_3_sel get y]
				set z3 [$atom_3_sel get z]
				set atom_4_id [lindex $ringListSingle $uniqueRingIndex 3]
				set atom_4_sel [atomselect top "index $atom_4_id" frame $n]
				set x4 [$atom_4_sel get x]
				set y4 [$atom_4_sel get y]
				set z4 [$atom_4_sel get z]
				set atom_5_id [lindex $ringListSingle $uniqueRingIndex 4]
				set atom_5_sel [atomselect top "index $atom_5_id" frame $n]
				set x5 [$atom_5_sel get x]
				set y5 [$atom_5_sel get y]
				set z5 [$atom_5_sel get z]
				set atom_6_id [lindex $ringListSingle $uniqueRingIndex 5]
				set atom_6_sel [atomselect top "index $atom_6_id" frame $n]
				set x6 [$atom_6_sel get x]
				set y6 [$atom_6_sel get y]
				set z6 [$atom_6_sel get z]
				
				# calculate the ring area
				set area_1 [area $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3]
				set area_2 [area $x5 $y5 $z5 $x4 $y4 $z4 $x3 $y3 $z3]
				set area_3 [area $x1 $y1 $z1 $x6 $y6 $z6 $x5 $y5 $z5]
				set area_4 [area $x1 $y1 $z1 $x3 $y3 $z3 $x5 $y5 $z5]
				set area_ring [expr {$area_1 + $area_2 + $area_3 + $area_4}]
				
				# calcualte the ring circularity
				set bond12 [measure bond [list $atom_1_id $atom_2_id] first $n last $n]
				set bond23 [measure bond [list $atom_2_id $atom_3_id] first $n last $n]
				set bond34 [measure bond [list $atom_3_id $atom_4_id] first $n last $n]
				set bond45 [measure bond [list $atom_4_id $atom_5_id] first $n last $n]
				set bond56 [measure bond [list $atom_5_id $atom_6_id] first $n last $n]
				set bond61 [measure bond [list $atom_6_id $atom_1_id] first $n last $n]
				set perimeter [expr {$bond12+$bond23+$bond34+$bond45+$bond56+$bond61}]
				set circularity [expr {4*3.141593*$area_ring/$perimeter/$perimeter}]
				set ringX [expr {($x1+$x2+$x3+$x4+$x5+$x6)/6.0}]
				set ringY [expr {($y1+$y2+$y3+$y4+$y5+$y6)/6.0}]
				set ringZ [expr {($z1+$z2+$z3+$z4+$z5+$z6)/6.0}]
				
				lappend uniqueRingAreaSingle $area_ring
				lappend uniqueRingCircularitySingle $circularity
				puts $ringAllFile [format "frameID %-8d index %-8d %-8d %-8d %-8d %-8d %-8d area  %-8.4f circularity  %-8.4f centerXYZ  %-6.2f %-6.2f %-6.2f" \
					  $frameID $atom_1_id $atom_2_id $atom_3_id $atom_4_id $atom_5_id $atom_6_id $area_ring $circularity $ringX $ringY $ringZ]
			}
		}
		
		set ringNumSingle [llength $uniqueRingListSingle]
		
		
		set ringAreaSumSingle 0
		foreach element $uniqueRingAreaSingle {
			set ringAreaSumSingle [expr {$ringAreaSumSingle + $element}]
		}
		
		set ringCircularitySumSingle 0
		foreach element $uniqueRingCircularitySingle {
			set ringCircularitySumSingle [expr {$ringCircularitySumSingle + $element}]
		}
		
		set avgRingAreaSingle [expr {$ringAreaSumSingle/$ringNumSingle}]
		set avgRingCircularitySingle [expr {$ringCircularitySumSingle/$ringNumSingle}]
		puts $ringAvgFile [format "%-8d %-8d %-8.4f %-8.4f" $frameID $ringNumSingle $avgRingAreaSingle $avgRingCircularitySingle]
		puts [format "%-8d %-8d %-8.4f %-8.4f" $frameID $ringNumSingle $avgRingAreaSingle $avgRingCircularitySingle]	
		close $ringAvgFile	
		close $ringAllFile
	} ; # ring mode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end

} ; # frame loop end

# for better visualization, wrap all atoms into the cell
# NOTE: the wraped bond connection may differ from the unwrapped bond connection, but just for visualization purpose
# pbc wrap -all
close $logfile

} ; # proc end


# subroutine for calculation of angle
proc angle {x1 y1 z1 x2 y2 z2 x3 y3 z3} {
	# Define the positions of the three points as vectors
	set vertex [list $x1 $y1 $z1]
	set point2 [list $x2 $y2 $z2]
	set point3 [list $x3 $y3 $z3]

	# Calculate vectors from the vertex to the other two points
	set vector1 [vecsub $point2 $vertex]
	set vector2 [vecsub $point3 $vertex]

	# Calculate the dot product of the two vectors
	set dot_product [vecdot $vector1 $vector2]

	# Calculate the magnitudes (lengths) of the vectors
	set length_1 [veclength $vector1]
	set length_2 [veclength $vector2]

	# Calculate the angle in radians using the dot product
	set cos_theta [expr {$dot_product / ($length_1 * $length_2)}]
	set angle_rad [expr {acos($cos_theta)}]

	# Convert the angle from radians to degrees
	set angle_degrees [expr {$angle_rad * 180.0 / 3.141593}]
	return $angle_degrees
}

# subroutine for calculation of ring area
proc area {x1 y1 z1 x2 y2 z2 x3 y3 z3} {
	# Define the positions of the three vertices (points) of the triangle as vectors
	set point1 [list $x1 $y1 $z1]
	set point2 [list $x2 $y2 $z2]
	set point3 [list $x3 $y3 $z3]

	# Calculate two vectors representing two sides of the triangle
	set vector1 [vecsub $point2 $point1]
	set vector2 [vecsub $point3 $point1]

	# Calculate the cross product of the two vectors
	set cross_product [veccross $vector1 $vector2]

	# Calculate the area of the triangle as half of the magnitude of the cross product
	set area [expr {[veclength $cross_product]/2.0}]
	return $area
}
