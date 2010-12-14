#############################################################
#  Ring polymers                                            #
#############################################################

puts "==================================================="
puts "=            Ring polymers                         ="
puts "==================================================="

puts "Program Information: \n[code_info]\n"

#############################################################
# File names and number of chains and monomers              #
#############################################################

# number of polymer chains
set M 200

# number of monomers per chain
set N 200

set box_l [expr pow(([expr $M*$N]/0.85),(1.0/3.0))]
setmd box_l $box_l $box_l $box_l

#############################################################
# Interactions                                              #
#############################################################

# repulsive Lennard Jones (Standard KG model)
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_shift   0.0
set lj1_off     0.0
set lj1_cut     1.12

# FENE
set fene_k     30.0
set fene_r      1.5

# angular potential
set angle_k     1.5

# define interactions
inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_off
inter 0   FENE          $fene_k  $fene_r
inter 1   angle         $angle_k
puts "Interactions:\n [inter]"


#############################################################
# Integration parameters                                    #
#############################################################

setmd skin 0.3
setmd time_step 0.01


#############################################################
# Read initial coordinates                                  #
#############################################################

set name "espresso_polymer_melt.start"
set inp [open "$name" r]
puts "Reading starting configurations ... "
while { [blockfile $inp read auto] != "eof" } {}
close $inp; puts -nonewline " done reading.\n"
flush stdout


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

puts "\n\n== System: $M chains of $N monomers =="
puts "Simulate $M polymer chains with $N monomers:"
puts "  box length: [setmd box_l]"
puts "  density: $density"
puts "  time step: [setmd time_step]"


#############################################################
# Langevin thermostat                                       #
#############################################################

set temperature 1.0
set gamma       1.0
thermostat langevin $temperature $gamma
#thermostat off
puts "Thermostat:\n [thermostat]"


#############################################################
# Perform integration                                       #
#############################################################

set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
set tmp_Temp [expr [analyze energy kin]/[setmd n_part]/1.5]
puts "  Analysis at t=[setmd time]: T=[setmd temp], P=$p1\n"
puts "\n\nEnergy==\n[analyze energy]"
puts "\n\nPressure tensor==\n[analyze stress_tensor]"; flush stdout

set msg [time { integrate 1000 } 1 ]
puts "  Elapsed time = $msg"

set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
set tmp_Temp [expr [analyze energy kin]/[setmd n_part]/1.5]
puts "  Analysis at t=[setmd time]: T=[setmd temp], P=$p1\n"
puts "\n\nEnergy==\n[analyze energy]"
puts "\n\nPressure tensor==\n[analyze stress_tensor]"; flush stdout

puts "Integration complete."; flush stdout; exit
