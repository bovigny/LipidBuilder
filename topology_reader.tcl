## Laboratory for Biomolecular Modeling, EPFL
## 2013 Version 1.0
## Thomas Lemmin

package provide topology_reader 1.0

namespace eval ::topology_reader:: {
    variable topology [list]
    variable impropertopology [list]
    variable atomtopology [list]
    variable chargetopology [list]
    variable bonds [list]
    variable ICs [list]
    variable masses [list]
    variable IC_linker [list]
    variable bonds_linker [list]
    variable limpropertopology [list]
}
proc ::topology_reader::init {} {
    variable topology 
    variable impropertopology 
    variable atomtopology 
    variable chargetopology
    variable bonds
    variable ICs
    variable masses 
    variable IC_linker
    variable bonds_linker 
    variable limpropertopology 
    set topology [list]
    set impropertopology [list]
    set atomtopology [list]
    set chargetopology [list]
    set bonds [list]
    set ICs [list]
    set masses [list]
    set IC_linker [list]
    set bonds_linker [list]
    set limpropertopology [list]

}

proc ::topology_reader::read_topology {file} {
        variable topology
        variable impropertopology
        variable atomtopology
        variable chargetopology
        variable masses
	array set masses_array $masses
        variable bonds 
	variable ICs 
        variable IC_linker
        variable bonds_linker
        variable limpropertopology
	set improperatoms [list]
        set improperlist [list] 
        set infile [open $file "r"]
        set resname "none"
        array set improper_array $impropertopology
	array set limproper_array $limpropertopology 
        while {[gets $infile line] >= 0} {
                if {[regexp {^RESI} $line]} {
                        set resname [lindex [noblanks [split $line ]] 1]
                }

                if {[regexp {^IMPR} $line]} {
                         set improperatoms [lreplace [noblanks [split $line ]] 0 0]
			 if {[info exists improper_array($resname)]} {
                                set b $improper_array($resname)
                                set b [concat $b $improperatoms]
                                set improper_array($resname) $b
                        } else {
                                set improper_array($resname) $improperatoms
                        }
                }
                
		if {[regexp {^ATOM} $line]} {
                        set atoms [lreplace [noblanks [split $line ]] 0 0]
			if {[info exists atom_array($resname)]} {
				set a $atom_array($resname)
				lappend a [lrange $atoms 0 2]
	                        set atom_array($resname) $a
			} else {
				set atom_array($resname) [list [lrange $atoms 0 2]]
			}
                        set charge_array($resname-[lindex $atoms 0]) [lindex $atoms 2]
                }
                if {[regexp {^MASS} $line]} {
                        set mass [lreplace [noblanks [split $line ]] 0 0]
                        set masses_array([lindex $mass 1]) [lrange $mass 2 end]
                }
                if {[regexp {^PRES} $line]} {
                        set resname "none"
                }
		if {[regexp {^IC} $line]} {
			set ic [lreplace [noblanks [split $line ]] 0 0]
			lassign $ic a1 a2 a3 a4 b a d aa bb
                        if {[info exists ic_array($resname)]} {
				set i $ic_array($resname)
				lappend i [list $a1 $a2 $a3 $a4]
				lappend i [list $b $a $d $aa $bb]
				set ic_array($resname) $i
			} else {
				set i [list [list $a1 $a2 $a3 $a4]]
				lappend i [list $b $a $d $aa $bb]
				set ic_array($resname) $i
			}
		}
		if {[regexp {^BOND} $line]|[regexp {^DOUBLE} $line]} {
			set bondatoms [lreplace [noblanks [split $line ]] 0 0]
			if {[info exists bond_array($resname)]} {
				set b $bond_array($resname)
				set b [concat $b $bondatoms]
				set bond_array($resname) $b
			} else {
				set bond_array($resname) $bondatoms
			}
				
		}
		if {[regexp {^LBOND} $line]} {
			set bondatoms [lreplace [noblanks [split $line ]] 0 0]
			if {[info exists lbond_array($resname)]} {
				set b $lbond_array($resname)
				set b [concat $b $bondatoms]
				set lbond_array($resname) $b
			} else {
				set lbond_array($resname) $bondatoms
			}
				
		}
		if {[regexp {^LIC} $line]} {
			set ic [lreplace [noblanks [split $line ]] 0 0]
			lassign $ic a1 a2 a3 a4 b a d aa bb
                        if {[info exists lic_array($resname)]} {
				set i $lic_array($resname)
				lappend i [list $a1 $a2 $a3 $a4]
				lappend i [list $b $a $d $aa $bb]
				set lic_array($resname) $i
			} else {
				set i [list [list $a1 $a2 $a3 $a4]]
				lappend i [list $b $a $d $aa $bb]
				set lic_array($resname) $i
			}
		}
                if {[regexp {^LIMPR} $line]} {
                         set improperatoms [lreplace [noblanks [split $line ]] 0 0]
			 if {[info exists limproper_array($resname)]} {
                                set b $limproper_array($resname)
                                set b [concat $b $improperatoms]
                                set limproper_array($resname) $b
                        } else {
                                set limproper_array($resname) $improperatoms
                        }
                }
        }
        close $infile
        set impropertopology [array get improper_array]
        set atomtopology [concat $atomtopology [array get atom_array]]
        set chargetopology [concat $chargetopology [array get charge_array]]
        set ICs [concat $ICs [array get ic_array]]
        set bonds [concat $bonds [array get bond_array]]
        set masses [concat $masses [array get masses_array]]
	set IC_linker [array get lic_array]
	set bonds_linker [array get lbond_array]
        set limpropertopology [array get limproper_array]
        set topology [list $atomtopology $chargetopology $masses $impropertopology $bonds $ICs $IC_linker $bonds_linker $limpropertopology]
}

proc ::topology_reader::read_linker {file} {
	set infile [open $file "r"]
	set resname "none"
        variable IC_linker
	while {[gets $infile line] >= 0} {
		if {[regexp {^RESI} $line]} {
			set resname [lindex [noblanks [split $line ]] 1]
		}
		if {[regexp {^IC} $line]} {
			set ic [lreplace [noblanks [split $line ]] 0 0]
                        lassign $ic a1 a2 a3 a4 b a d aa bb
                        if {[info exists ic_array($resname)]} {
                                set ic_array($resname) [lappend $bond_array($resname) [list $a1 $a2 $a3 $a4]]
                                set ic_array($resname) [lappend $bond_array($resname) [list $b $a $d $aa $bb]]
                        } else {
                                set ic_array($resname) [list $a1 $a2 $a3 $a4]
                                set ic_array($resname) [list $b $a $d $aa $bb]
                        }
                }
	}
	close $infile
	set IC_linker [array get ic_array]
}

proc ::topology_reader::noblanks {mylist} {
  set newlist [list]
  foreach elem $mylist {
    if {$elem != ""} {
      lappend newlist $elem
    }
  }

  return $newlist
}
