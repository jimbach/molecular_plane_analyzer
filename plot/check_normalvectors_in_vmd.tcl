proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.15
    graphics $mol cone $middle $end radius 0.3
}





proc checkallvectors {} {
  set numberofatoms [molinfo top get numatoms]		;# all atom indices
  set allatomlist [[atomselect top "all"] get type] 	;# all atom types
  set allcoordlist [[atomselect top "all"] get {x y z}]	;# all atom coords
  
  for {set i 0} {$i < $numberofatoms} {incr i} {
	set activesite [lindex $allcoordlist $i]
	set activetype [lindex $allatomlist $i]
	if { $activetype == "C" } then {
		set fromatom $activesite
	}
	
	if { $activetype == "O" } then {
		draw color red
		draw arrow $fromatom  $activesite
	}
	
	if { $activetype == "N" } then {
		draw color blue
		draw arrow $fromatom $activesite
	}
  }
}