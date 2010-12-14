#############################################################
#  Lennard-Jones NVE                                        #
#############################################################

puts "==================================================="
puts "=             Lennard-Jones NVE                   ="
puts "==================================================="

puts "Program Information: \n[code_info]\n"

#############################################################
# Read initial coordinates                                  #
#############################################################

set name "espresso_lennard_jones.start"
set inp [open "$name" r]
puts "Reading starting configurations ... "
while { [blockfile $inp read auto] != "eof" } {}
close $inp; puts -nonewline " done reading.\n"
flush stdout


#############################################################
# Interactions                                              #
#############################################################

# repulsive Lennard Jones (Standard KG model)
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_shift   0.0
set lj1_off     0.0
set lj1_cut     2.5

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_off
puts "Interactions:\n [inter]"


#############################################################
# Integration parameters                                    #
#############################################################

set steps 1000
setmd skin 0.3
setmd time_step 0.005


#############################################################
# Other parameters                                          #
#############################################################

set tcl_precision 5
set random_seeds { }


#############################################################
# Random number generator setup                             #
#############################################################

if { [llength $random_seeds] > 0 } {
   eval t_random seed $random_seeds
}


#############################################################
# Write out values                                          #
#############################################################

puts "\n\n== System:  =="
puts "Simulate LJ fluid:"
puts "  box length: [setmd box_l]"
puts "  time step: [setmd time_step]"
puts "  density: $density"


#############################################################
# Langevin thermostat                                       #
#############################################################

set temperature 1.0
set gamma       1.0
#thermostat langevin $temperature $gamma
thermostat off
puts "Thermostat:\n [thermostat]"


#############################################################
# Perform integration                                       #
#############################################################

set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
set tmp_Temp [expr [analyze energy kin]/[setmd n_part]/1.5]
puts "  Analysis at t=[setmd time]: T=[setmd temp], P=$p1\n"
puts "[analyze energy]"
puts "[analyze stress_tensor]"; flush stdout

set start_time [clock seconds]
set msg [time { integrate $steps } 1 ]
set end_time [clock seconds]
puts "TCL time command = $msg"
puts "TCL clock = [expr $end_time - $start_time]"

set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
set tmp_Temp [expr [analyze energy kin]/[setmd n_part]/1.5]
puts "\nAnalysis at t=[setmd time]: T=[setmd temp], P=$p1\n"; flush stdout
puts "[analyze energy]"
puts "[analyze stress_tensor]"; flush stdout

puts "Integration complete."; flush stdout; exit
