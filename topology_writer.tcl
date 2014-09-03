## Laboratory for Biomolecular Modeling, EPFL
## 2013 Version 1.2
## Thomas Lemmin

package provide topology_writer 1.2

namespace eval ::topology_writer:: {
    variable atoms [list]
    variable bonds [list]
    variable ICs [list]
    variable masses [list]
    variable impropers [list]
   
}

proc ::topology_writer::init {} {
	variable atoms
	set atoms [list]
	variable bonds 
	set bonds [list]
	variable ICs 
	set ICs	[list]
	variable masses 
	set masses [list]
	variable impropers 
	set impropers [list]
}

proc ::topology_writer::write_topology {file resname charge} {
	variable atoms
	variable bonds
	variable impropers 
	variable ICs
	variable masses
	set outfile [open $file "w"]
	#HEADER
	puts $outfile {*  \\\\\\\ CHARMM36 All-Hydrogen Lipid Topology File ///////}
	puts $outfile {*  All comments and questions should be submitted to the}
	puts $outfile {*  parameter forum at the CHARMM website: www.charmm.org}
	puts $outfile {*}
	puts $outfile ""
	puts $outfile "36  1"
	puts $outfile ""
	puts $outfile "! Created with LipidBuilder 2.1"
	puts $outfile ""
	puts $outfile ""
	#MASS
	set i 1
	foreach {k m} $masses {
		puts $outfile "MASS $i $k\t$m"
		incr i
	}
	puts $outfile ""
	puts $outfile ""
	puts $outfile ""
	puts $outfile "RESI $resname  $charge"
	#ATOM
	foreach a $atoms {
		puts $outfile "ATOM $a"
	}
	puts $outfile ""
	#BONDS
	set counter 1
	puts -nonewline $outfile "BOND "
	foreach {b1 b2} $bonds {
		if {[expr $counter%5]==0} {
			puts -nonewline $outfile "\nBOND $b1 $b2 "
		} else {
			puts -nonewline $outfile "$b1 $b2 "
		}
		incr counter
	}
	
	puts $outfile ""
	puts $outfile ""
	foreach {a1 a2 a3 a4} $impropers {
		puts $outfile "IMPR $a1 $a2 $a3 $a4"
	}
	puts $outfile ""
	puts $outfile ""
	puts $outfile "!Internal coordinates"
	#IC
	foreach {a ic} $ICs {
		puts $outfile "IC $a $ic"
	}
	puts $outfile ""
	puts $outfile ""
	puts $outfile "AUTO ANGLES DIHE"
	puts $outfile ""

	close $outfile
}

proc ::topology_writer::concat_top {hatoms hbonds hICs} {
	variable ICs
	variable bonds
	variable atoms
	set ICs [concat  $ICs $hICs]
	set bonds [concat $bonds $hbonds]
	set atoms [concat $atoms $hatoms]
}

proc ::topology_writer::set_masses {topo_masses} {
	variable atoms
	variable masses
	array set tmasses $topo_masses
	set aunique [lsort -unique -index 1 $atoms]
	foreach a $aunique {
		lassign $a n t c
		set masses_array($t) $tmasses($t)
	}
	
	set masses [array get masses_array]
}

proc ::topology_writer::set_connectivity_IC {lIC IC_linker} {
    variable ICs
	foreach {a ic} $IC_linker {
		set aa $a
		set i 1
		foreach c $lIC {
			set c [split $c ""] 
			#lIC are sorted: C H H C H H ... 
			if {[lindex $c 0] eq "C"} {
				incr i
				set l "L[lindex $c 1]$i"
			} else {
				set l "P${i}[lindex $c end]"
			}
			set l [join $l ""]
			set c [join $c ""]
			set aa [string map [list $l $c] $aa]
		}
		#check if linker has been changed
		if {![string match *L* $aa] && ![string match *P* $aa]} {
			lappend ICs $aa
			lappend ICs $ic
		} else {
			#puts "Invalid Linker $aa"
		}
	}
}

proc ::topology_writer::set_connectivity_bonds {lIC lbonds} {
	variable bonds
	foreach {b1 b2} $lbonds {
		set b [list $b1 $b2]
		foreach c $lIC {
			set b [replace_linker $b $c]
		}
		if {![string match *L* $b] && ![string match *P* $b]} {
			set bonds [concat $bonds $b]
		} else {
			#puts "Invalid bond "
		}
	}
}

proc ::topology_writer::set_improper {lIC limpropers} {
	variable impropers
	foreach {a1 a2 a3 a4} $limpropers {
		set a [list $a1 $a2 $a3 $a4]
                set aa $a
                foreach c $lIC {
					set aa [replace_linker $aa $c]
                }
                #check if linker has been changed and no dummy atoms (L or P) remain
                if {[string compare [join $a] [join $aa]] && ![string match *L* $aa] && ![string match *P* $aa]} {
                        lappend impropers $aa
                } else {
                        #puts "Invalid improper $aa"
                }
        }
}

#Replace L is used to define carbon atom (C) in linker and P for hydrogen (H)
#=>replaces L* with C* and P* with H*
proc ::topology_writer::replace_linker {link aa} {
			set linkaa [string map {H P} [string map {C L} $aa]]
			return [string map [list $linkaa $aa] $link]
}
