package require readcharmmtop

proc ::Topowrite::write_charmm_topology { top file {header {}} {resi {}} } {
   set fd [open $file "w"]
   #HEADER
   if {![llength $header} {
      set header [list {*>>>>>> CHARMM topology file generated by LipidBuilder <<<<<<}\
                       {36 1}]
   }
   foreach h $header {
     puts $fd $h
   }
   puts $fd ""

   #format all residue
   if {![llength $resi]} {
       set resi [::Toporead::topology_get names $top]
   }
   set topcontents [list]
   set masscontents [list]
   set masses [::Toporead::topology_get types $top]
   #RES
   foreach r $resi {
      set name [::Toporead::topology_get_resid $top $r name]
      set charge [::Toporead::topology_get_resid $top $r charge]
      set head [::Toporead::topology_get_resid $top $r head]
      lappend topcontents [format "%4s %4s    % 5.2f\n" $head $name $charge]
      #GROUP 
      foreach group [::Toporead::topology_get_resid $top $r grouplist] {
         lappend topcontents "GROUP\n"
         foreach atom $group {
            lassign $atom name type charge
            lappend topocontents [format "ATOM %4s %4s % .5f\n" $name $type $charge]
            lappend masscontents [lsearch -inline -index 0  $masses $type]
          }
       }
       lappend topocontents "\n"
       #DELETE for PRES
       foreach d [::Toporead::topology_get_resid $top $r delete] {
          lappend topocontents [format "DELETE ATOM $4s\n" $d]
       }
       lappend topocontents "\n"
       #BOND
       lappend topocontents "BOND "
       set nbonds 0
       foreach bond [::Toporead::topology_get_resid $top $r bonds] {
           lappend topocontents "[join $bond " "] "
           if {[expr $nbonds % 4] == 0} {
              lappend topcontents "\n"
              lappend topcontents "BOND "
            }
            incr nbonds
       }
       lappend topocontents "\n"
       #IMPR
       foreach imp [::Toporead::topology_get_resid $top $r impropers] {
           lappend topocontents "IMPR [join $imp " "]"
        }
       lappend topocontents "\n"
       #IC
       foreach ic [::Toporead::topology_get_resid $top $r ics] {
           lassign $ic atom1 atom2 atom3 atom4 bond1 angle1 dihedral1 angle2 bond2
           lappend topocontents [format "IC %4s %4s %4s %4s  % .5f % .5f % .5f % .5f % .5f\n"\ 
                                 $atom1 $atom2 $atom3 $atom4 $bond1 $angle1 $dihedral1 $angle2 $bond2]
        }
        lappend topcontents "\n\n"
   }
   #MASS
   set idx 1
   foreach mass $masscontents {
        lassign $mass type mass elem comment
        puts $fd [format "MASS %5i %-4s  %8.5f %2s !%s" $idx $type $mass $elem $comment]
   }
   puts $fd ""
   puts $fd "AUTO ANGLES DIHE"
   puts $fd ""  
   }
   #write all the residues/patch 
   foreach line $topocontents {
      puts -nonewline $fd $line
   }
   puts $fd "\nEND\n\n"
   close $fd
}