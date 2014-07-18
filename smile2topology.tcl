# Assigns atom type and internal coordinates of a SMILE sequence, in order to build the topology file and 3D model.
#
#
# Thomas Lemmin
#	Laboratory for Biomolecular Modeling, EPFL, Switzerland
#	2013
#

package provide smile2topology 1.0
package require smiles_parser 1.0

namespace eval ::smile2topology:: {
	
	set numElkey 6

		
	#key in ICkey
	set ielement 0
	set ilinker 1
	set iconnectivity 2
	set iring 3
	set ibranch 4
	set ibond 5
	
	#value in atoms
	set iname 0
	set itype 1
	set icharge 2
	set ihydrogen 4
	#resi list
	set rname 0
	set rhatom_name 1
	set rkey {2 0}
	set ratom {2 1}
	set rsterep {2 2}
	#variable for topology writer
    variable atoms [list]
    variable bonds [list]
    variable ICs [list]
    
    #private variable
    variable heavy_atom_key [list]
    variable heavy_atom_name [list]
    variable ICparameters [list]
    variable ICkey [list]
    variable ICresi [list]    
    variable ICbond [list]    
    variable ICangle [list]
    variable ICdihedral [list]

}

#initialize all the variables
proc ::smile2topology::init {} {
	variable atoms
    variable bonds
    variable ICs
    #private variable
    variable heavy_atom_key
    variable heavy_atom_name
    variable hydrogens
    variable ICparameters
    variable ICresi
    variable ICbond
    variable ICangle
    variable ICdihedral

	set atoms [list]
    set bonds [list]
    set ICs [list]
    #private variable
	set heavy_atom_key [list]
	set heavy_atom_name [list]
    set ICparameters [list]
    set ICkey [list]
    set ICresi [list]
    set ICbond [list]
    set ICangle [list]
    set ICdihedral [list]
}

proc ::smile2topology::initsmile {} {
	variable atoms
    variable bonds
    variable ICs
    #private variable
    variable heavy_atom_key
    variable heavy_atom_name
    variable hydrogens
    set atoms [list]
    set bonds [list]
    set ICs [list]
    #private variable
	set heavy_atom_key [list]
	set heavy_atom_name [list]
}

proc ::smile2topology::noblanks {mylist} {
  set newlist [list]
  foreach elem $mylist {
    if {$elem != ""} {
      lappend newlist $elem
    }
  }

  return $newlist
}

#Reads IC parameters
#	RESI [defines atomtype, charge and hydrogens]
#	 	name of resi is used to identify BOND, ANGLE, DIHEDRAL, STEREO
#		KEY: smile representation of RESI:
#					ELEMENT LINKER CONNECTIVITY RING BRANCH BOND
#					ELEMENT: chemical element [C, N, O, ...]
#					LINKER: first atom in chain [0; 1]
#					CONNECTIVITY: number of heavy atoms bound to atom
#					RING: atom in a ring [0; 1]
#					BOND: bond type (single, double, triple bond): [1; 2; 3] (e.g C-C-C = 1:1)
#		ATOM 
#			atoms composing a RESI
#			NAME CHARGE TYPE
#		HATOM
#			heavy atom of RESI
#			NAME CHARGE TYPE
#		STEREO [defines the stereochemistry between four atoms in a RESI]
#			ATOM1 ATOM2 ATOM3 ANGLE (e.g. CH2 CH2 CH2 CH3 120)
#
#	BOND [defines the theoretical distance between two atoms]
#		RESI:ATOM1 RESI:ATOM2 BOND (e.g. CH3:C CH2:C 1.5)
#	ANGLE [defines the theoretical angle between three atoms]
#		RESI:ATOM1 RESI:ATOM2 RESI:ATOM3 RESI:ANGLE (e.g. CH2:C CH1:C CH1:H 120)
#	DIHEDRAL [defines the theoretical angle between four atoms]
#		RESI:ATOM1 RESI:ATOM2 RESI:ATOM3 RESI:ATOM4 ANGLE (e.g. CH2:C CH2:C CH2:C CH3:C 180)
#
#	* are used to indiquate
#	! are used for comments
#

proc ::smile2topology::read_ICparameters {file} {
	variable ICparameters
	variable ICkey
	variable ICresi
	variable ICbond 
	variable ICangle
	variable ICdihedral
	
	set type "none"
	set resname "none"
	set keys [list]
	set atoms [list]
	set infile [open $file "r"]
        while {[gets $infile line] >= 0} {
                switch -regexp  $line {
			{^\s*RESI} {
			set type "resi"
			}
			{^\s*ISER} {
			#save residue
			set ICresi_array($resname) [list [list "key" $keys] [list "atom" $atoms] 	[array get ICstereo_array]]
			# reinitialize variable and lists
			set type "none"
			set resname "none"
			set keys [list]
			set atoms [list]
			if { [array exists ICstereo_array] } {
				array unset ICstereo_array
			}
			continue}
			{^\s*KEY} {
			set type "key"
			}
			{^\s*ATOM} {
			set type "atom"
			}
			{^\s*HATOM} {
			set type "hatom"
			}
			{^\s*STEREO} {
			set type "stereo"
			}
			{^\s*BOND} {
			set type "bond"
			continue}
			{^\s*ANGLE} {
			set type "angle"
			continue}
			{^\s*DIHEDRAL} {
			set type "dihedral"
			continue}
		}
		if {[string compare $type "none"] != 0} {
			set linearray [split [lindex [split $line !] 0 ]]
        	set linearray [noblanks $linearray]
             	if {[llength $linearray] > 0} {
					switch $type {
						"resi" {
							set resname [lindex $linearray 1]
						}
						"key" {
							set linearray [lrange $linearray 1 end]
							if {[llength $linearray] == $::smile2topology::numElkey} {
								lappend keys [join $linearray =]
								set ICkey_array([lindex $keys end]) $resname
						} else {
								puts "WARNING invalid KEY in resi $resname: $linearray"
							}
						}
						"atom" {
							set linearray [lrange $linearray 1 end]
							if {[llength $linearray] == 3} {
								lappend atoms [concat "hydrogen" $linearray]
							} else {
								puts "WARNING invalid ATOM in resi $resname: $linearray"
							}
						}
						"hatom" {
							set linearray [lrange $linearray 1 end]
							if {[llength $linearray] == 3} {
								lappend atoms [concat "heavy" $linearray]
							} else {
								puts "WARNING invalid ATOM in resi $resname: $linearray"
							}
						}
						"stereo" {
							set linearray [lrange $linearray 1 end]
							if {[llength $linearray] >= 5} {
								set b [lassign $linearray a1 a2 a3 a4]
								set k $a1=$a2=$a3=$a4
								if {[info exists ICstereo_array($k)]} {
									set aa $ICstereo_array($k)
									set b [concat $aa $b] 
								}
								set ICstereo_array($k) $b
							} else {
								
								puts "WARNING invalid STEREO in resi $resname: $linearray"
							}
						}
						"bond" {
							if {[llength $linearray] == 3} {
								lassign $linearray a1 a2 b
								set ICbond_array($a1=$a2) $b
							} else {
								puts "WARNING invalid bond: $linearray"
							}
						}
						"angle" {
							 if {[llength $linearray] == 4} {
													lassign $linearray a1 a2 a3 b
													set ICangle_array($a1=$a2=$a3) $b
											} else {
													puts "WARNING invalid angle: $linearray"
											}
						}
						"dihedral" {
							if {[llength $linearray] >= 5} {
								set b [lassign $linearray a1 a2 a3 a4]
								set ICdihedral_array($a1=$a2=$a3=$a4) $b
							} else {
								 puts "WARNING invalid dihedral: $linearray"
							}
						}
        		}
  			}
		}
	}
	close $infile
	#Convert arrays to lists
	set ICresi [array get ICresi_array]
	set ICkey [array get ICkey_array]
	set ICbond [array get ICbond_array]
	set ICangle [array get ICangle_array]
	set ICdihedral [array get ICdihedral_array]
	set ICparameters [list $ICresi $ICbond $ICangle $ICdihedral]
}	


#########################################################################
#					assign connectivity of atoms						#
#########################################################################



#parse the SMILE flatten structure and assigns the atom key
proc ::smile2topology::assign_atom_keys {structure} {
	variable heavy_atom_key
	
	lassign $structure atom bond chains
	
	foreach a $atom {
		set aelement [lindex $a $::smiles_parser::atom_symbol];
		set alinker 0;
		if {[lindex $a $::smiles_parser::atom_id] == 0} {
			set alinker 1
		}
		set bondlist [lindex $a $::smiles_parser::atom_linkedby]
		set aconnectivity [llength $bondlist]
		set abonds [list]
		set aring 0
		set abranch 0
		if {[llength [lindex $a $::smiles_parser::atom_branches]]> 0} {
			set abranch 1
		}
		
		foreach b $bondlist {
			lappend abonds [lindex [lindex $bond $b] $::smiles_parser::bond_count]
		}
		
		set key [list $aelement $alinker $aconnectivity $aring $abranch $abonds]
		lappend heavy_atom_key $key
	}
	correct_ringatoms $structure
	correct_isomer $structure
}		


proc ::smile2topology::detect_ringatoms {structure} {
	set index [list]
	lassign $structure atom bond chain
	
	#find all ring bond
	foreach c $chain {
		set ringatoms [lsearch -all [lindex $c $::smiles_parser::chain_list] *ringbond*]
		foreach r $ringatoms {
			#start from first atom of ring
			set sourceatom [lindex [lindex $c $::smiles_parser::chain_list] $r]
			set ringbond [lindex $sourceatom $::smiles_parser::atom_ringbonds]
			if {[llength $ringbond]>0} { 
				set bond_in_chain [list]
				#save all bonds in chain
				foreach e [lindex $c $::smiles_parser::chain_list] {
					if {[lindex $e 0] == "bond" } {
						lappend bond_in_chain [lindex $e $::smiles_parser::bond_id]
					}
				}
				#extract ringbond
				set ringbond [lindex [lindex [join $ringbond] $::smiles_parser::ringbond_bond] $::smiles_parser::bond_linking]
				set sourceatom [lindex $ringbond 0]
				lappend index $sourceatom
				set targetatom [lindex $ringbond 1]
				set currentatom $sourceatom
				
				#find path connecting source to target
				set count 0
				while { $currentatom != $targetatom } {			
					set bondlist [lindex [lindex $atom $currentatom] $::smiles_parser::atom_linkedby]
					#check if bond belongs to chain
					foreach currentbond $bondlist {
						if {[lsearch $bond_in_chain $currentbond] !=-1 } {
							set currentatom [lindex [lindex [lindex $bond $currentbond] $::smiles_parser::bond_linking] 0]
							lappend index $currentatom
							break
						}	
					}
					incr count
					if {$count>5} {
						puts "Oooops, I got lost trying to connect atom $sourceatom to atom $targetatom"
						break
					}
				}
			}
		}
	}
	
	return $index
}

proc ::smile2topology::correct_ringatoms {structure} {
	variable heavy_atom_key
	set index [detect_ringatoms $structure]
	foreach i $index {
		set key [lindex $heavy_atom_key $i]
		set key [lreplace $key $::smile2topology::iring $::smile2topology::iring 1]
		set heavy_atom_key [lreplace $heavy_atom_key $i $i $key]
	}
}

proc ::smile2topology::correct_isomer {structure} {
	variable heavy_atom_key
	lassign $structure atom bond chain
	foreach c $chain {
		set isomer [lsort -integer [concat [lsearch -all [lindex $c $::smiles_parser::chain_list] */*] [lsearch -all [lindex $c $::smiles_parser::chain_list] *\\\\*]]]
		if { [llength $isomer]>1 } {
			set c1 [lindex $isomer 0]
			set isomer [lreplace $isomer 0 0]
			foreach c2 $isomer {
				set b1 [lindex [lindex $c $::smiles_parser::chain_list] $c1]
				set b2 [lindex [lindex $c $::smiles_parser::chain_list] $c2]
				set a1 [lindex [lindex $b1 $::smiles_parser::bond_linking] 1]
				set a2 [lindex [lindex $b2 $::smiles_parser::bond_linking] 0]
				set abond [list $a1 $a2]
				if {[lsearch $bond "*$abond*"] !=-1} {
				#cis bond
					if {[lindex $b1 $::smiles_parser::bond_direction] != [lindex $b2 $::smiles_parser::bond_direction]} {
						set isomer1 [lindex [lindex [lindex $heavy_atom_key [lindex [lindex $b1 $::smiles_parser::bond_linking] 1]] $::smile2topology::ibond] 1]
						set i [lindex [lindex $b2 $::smiles_parser::bond_linking] 0]
						set key2 [lindex $heavy_atom_key $i]
						set isomer2 [lindex $key2 $::smile2topology::ibond]
						set isomer2 [lreplace $isomer2 1 1 [expr -1*($isomer1/abs($isomer1))*[lindex $isomer2 1]]] 
						set key2 [lreplace $key2 $::smile2topology::ibond $::smile2topology::ibond $isomer2]
						set heavy_atom_key [lreplace $heavy_atom_key $i $i $key2]
					}
				}
				set c1 $c2
			}
		}
	}
}


#########################################################################
#				proc for assigning topology parameters					#
#########################################################################

proc ::smile2topology::guess_quadruplet {structure} {
	lassign $structure atom bond chain
	set quadruplet [list]
    foreach b $bond {
    	lassign [lsort -integer [lindex $b $::smiles_parser::bond_linking]] a1 a2
    	set bond1 [lindex $atom [list $a1 $::smiles_parser::atom_linkedby]]
    	set bond2 [lindex $atom [list $a2 $::smiles_parser::atom_linkedby]]
    	if {[llength $bond1] > 1 && [llength $bond2] > 1} {
			foreach b1 $bond1 {
				set o1 [lsearch -not -inline [lindex $bond [list $b1 $::smiles_parser::bond_linking]] $a1]
				foreach b2 $bond2 {
					set o2 [lsearch -not -inline [lindex $bond [list $b2 $::smiles_parser::bond_linking]] $a2]
					if {($o1 == $a1) || ($o2 == $a1) || ($o1 == $a2) || ($o2 == $a2) ||($o2 == $o1) || ($o1 > $a1) || ($o2 < $a2) } {
                    	continue
                	}
					lappend quadruplet [list $o1 $a1 $a2 $o2]    		
				}
			}
    	}
    }
   return $quadruplet
}


proc ::smile2topology::assign_atoms {prefix_C suffix_H} {
	variable heavy_atom_key
	variable heavy_atom_name
	variable atoms
    variable ICresi
    variable ICkey
    #initialise list of atom
    set atoms [list]
    set heavy_atom_name [list]
    array set ICkey_array $ICkey
    array set ICresi_array $ICresi
    set index 1
	foreach key $heavy_atom_key {
		incr index
		set bonds [lindex $key $::smile2topology::ibond]
		set key [join [lreplace $key $::smile2topology::ibond $::smile2topology::ibond [join $bonds ":"]] "="]
		if {[info exists ICkey_array($key)]} {
			set resname $ICkey_array($key)
			lassign $ICresi_array($resname) rkey ratom rstereo 
			set ratom [lindex $ratom 1]
			set nameC "${prefix_C}${index}"
			lassign [lsearch -inline $ratom "*heavy*"] tag namec typeC chargeC
			set  hydrogens [lsearch -all -inline $ratom "*hydrogen*"]
			set hydrogenC [list]
			foreach h $hydrogens sh $suffix_H {
				if {$h=="" || $sh==""} {
					break
				}
				lassign $h t n t c
				lappend hydrogenC [list  $n "H${index}$sh" $t $c]
			}			
			lappend atoms [concat [list $namec $nameC $typeC $chargeC] [list $hydrogenC]]
			lappend heavy_atom_name "$resname:$namec:$nameC"
		} else {
			puts "ERROR, unknown atom $key"
			return -1
		}
	}
	return 1
}

proc ::smile2topology::assign_bonds {structure} {
	lassign $structure atom bond chain
	variable atoms
	variable bonds
	set nC [list]
	foreach b $bond {
		lassign [lindex $b $::smiles_parser::bond_linking] a1 a2
		lappend bonds [list [lindex $atoms [list $a1 1]] [lindex $atoms [list $a2 1]]]
		foreach a [list $a1 $a2] {
			set nameC [lindex $atoms [list $a 1]]
			set H [lindex $atoms [list $a $::smile2topology::ihydrogen]]
			foreach h $H {
				lappend bonds [list $nameC [lindex $h 1]]
			}
		}
	}
	set bonds [lsort -dictionary -unique $bonds]
}

proc ::smile2topology::get_linker {} {
	variable atoms
	variable bonds
	set atom [lindex $atoms {1 1}]
	set latom [lsort -unique [join [lsearch -inline -all -regexp $bonds "$atom C|C* $atom"] " "]]
	
	set atom_topo [list]
	foreach la $latom {
		set a [lsearch -inline $atoms *$la*]
		lappend atom_topo [lindex $a 1 ]
			foreach h [lindex $a $::smile2topology::ihydrogen] {
				lappend atom_topo [lindex $h 1]
			}
		}
		return $atom_topo
}

proc ::smile2topology::get_atoms {} {
		variable atoms
		set atom_topo [list]
		foreach a $atoms {
			lappend atom_topo [lrange $a 1 3]
			foreach h [lindex $a $::smile2topology::ihydrogen] {
				lappend atom_topo [lrange $h 1 3]
			}
		}
		return $atom_topo
}

proc ::smile2topology::get_bonds {} {
	variable bonds
	return [join $bonds " "]
}

proc ::smile2topology::set_ICs {structure} {
	variable ICs
	variable atoms
	variable ICresi
	variable ihydrogen
	variable heavy_atom_name
	variable heavy_atom_key
	set quadruplet [guess_quadruplet $structure]
	set count 0
	array set ICresi_array $ICresi
	foreach an $heavy_atom_name {
		set quadatom [lsearch -all -inline -regexp $quadruplet "^\\d+ \\d+ $count \\d+$"]
		set lquad [llength $quadatom]
		switch $lquad {
			#Terminal atom?
			0 {
				set quadatom [lsearch -inline -exact -regexp $quadruplet "^\\d+ \\d+ \\d+ $count\$"]
				if {[llength $quadatom] > 0 } {
					set key [list]
					lassign $quadatom a1 a1 a2 a3
					set quadatom [list $a1 $a2 $a3 $a3]
					foreach aa $quadatom {
						lappend key [join [lrange [split [lindex $heavy_atom_name $aa] ":" ] 0 1 ] ":"]
					}
					set names [list [lindex $atoms [list $a1 1]] [lindex $atoms [list $a2 1]] [lindex $atoms [list $a3 1]] [lindex $atoms [list $a3 1]]]
					#names of heavy atoms
					set stereo [lindex $ICresi_array([lindex [split [lindex $heavy_atom_name $a3] ":"] 0]) 2]
					set hydrogens []
					foreach h [lindex $atoms $count $ihydrogen] {
						lappend hydrogens [lrange $h 0 1]
					}
					#look stereo
					set_stereo $key $names $hydrogens $stereo
				}
			} default {
				#Standard atom
				set key [list]
				set names [list]
				set iso [list]
				set stereo [lindex $ICresi_array([lindex [split [lindex $heavy_atom_name [lindex $quadatom {0 2}]] ":"] 0]) 2]
				set hydrogens []
				foreach h [lindex $atoms $count $ihydrogen] {
					lappend hydrogens [lrange $h 0 1]
				}	
				foreach quada $quadatom {
					set quada [join $quada]
					lassign $quada a1 a2 a3 a4
					lappend names [list [lindex $atoms [list $a1 1]] [lindex $atoms [list $a2 1]] [lindex $atoms [list $a3 1]] [lindex $atoms [list $a4 1]]]
					set k [list]
					foreach aa $quada {
						lappend k  [join [lrange [split [lindex $heavy_atom_name $aa] ":" ] 0 1 ] ":"]
					}
					lappend key $k
					set iiso [lindex $heavy_atom_key [list $a3 5 end]]
					lappend iso [expr $iiso/abs($iiso)]
					#look stereo
					set_stereo [lindex $key end] [lindex $names end] $hydrogens $stereo
				}
				#look up IC
				set_topology $key $names $iso

			} 
		}
		incr count
	}
}


proc ::smile2topology::set_topology {key names iso} {
	variable ICbond
	variable ICangle
	variable ICdihedral
	variable ICs
	set d1 0
	foreach k $key n $names {
		set di  [expr [lindex [lookup_topology [join $k "="] $ICdihedral] $d1] - 90*(1-[lindex $iso $d1]) ]
		set bd1 [lookup_topology [join [lrange $k 0 1] "="] $ICbond]
		set ae1 [lookup_topology [join [lrange $k 0 2] "="] $ICangle]
		set ae2 [lookup_topology [join [lrange $k 1 3] "="] $ICangle]
		set bd2 [lookup_topology [join [lrange $k 2 3] "="] $ICbond]
		lappend ICs $n
		lappend ICs [list $bd1 $ae1 $di $ae2 $bd2]
		set d1 end
	}

}

proc ::smile2topology::set_stereo {key names hydrogens table} {
	variable ICs
	array set table_array $table
	set lkey [split $key "="]
	lassign [split [lindex $key 2] ":"] resi atom
	set stereocenter "${resi}:${atom}"
	set centernames [lindex $names 2]
	foreach hh $hydrogens {
		lassign $hh h hn
		lappend stereocenter "$resi:${h}"
		lappend centernames $hn
	}
	set lstereo [list]
	foreach  k [array names table_array] {
		set nstereo [list]
		foreach sk  [split $k "="] {
			set atomi [expr 2 + [llength [lsearch -all [split $sk ""] "+"]]]
			if {$atomi == 2} {
				set atomi [expr 2 - [llength [lsearch -all [split $sk ""] "-"]]]				
			}
			lassign [split $sk ":"] kresi katom
			if {$atomi == 2} {
				set cn 0
				foreach sc $stereocenter {
					if {$sk eq $sc || $sk eq "&$sc"} {
						set improper ""
						if {[string range $sk 0 0] eq "&"} {
							set improper "*"
						}
						lappend nstereo "$improper[lindex $centernames $cn]"
						break
					}
					incr cn
				}
			} else {
				lassign [split [lindex $lkey [list 0 $atomi]] ":"] resi atom
				if {($kresi eq "*" || $kresi eq $resi) && ([lindex [split $katom ""] 0] eq [lindex [split $atom ""] 0] || $katom eq $atom)} {
					lappend nstereo [lindex $names $atomi]
				} else {
					set nstereo "none"
					break
				}
			}
		}
		if {$nstereo != "none" && [llength $nstereo] == 4} {
			lappend ICs [join $nstereo " "] 
			lappend ICs $table_array($k)
		} else {
			puts "$k"
		}
	}
	return $lstereo			
}

proc ::smile2topology::lookup_topology {key table} {
	variable ICparameters
	array set table_array $table
	set table_keys [array names table_array]
	#check if key exists
	if {[info exists table_array($key)]} {
			return $table_array($key)
	}
	#search for closest key
	#build regular expression
	set key_rEx "(\\*|[regsub -all ":(\\w)" [regsub -all "=" $key {|\\*)=(\\*|}] {):\1(}]|\\*)"
	#search for all valid keys
	set valid_keys [lsearch -all -inline -regexp $table_keys $key_rEx]
	#get the closest one to original key (i.e. with the least *)
	if {[llength $valid_keys]>0} {
		set lkey [list]
		foreach vk $valid_keys {
			lappend lkey [regexp -all {\*} $vk]
		}
		set ikey [lsort -indices $lkey]
		set thekey [lindex $valid_keys [lindex $ikey 0]]
		return $table_array($thekey)
	} else {
		puts "WARNING: The parameter for [split $key =] was not found in topology"
		return ""
	}
	
}
