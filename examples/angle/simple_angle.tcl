setmd box_l 10.0 10.0 10.0
setmd skin 0.3
setmd time_step 0.001
thermostat off

inter 0 0 lennard-jones 1.0 1.0 2.5 0.0 0.0
inter 0   fene 30.0 1.5
inter 1   angle 1.5 3.1415926

part 0 pos 5.0 5.0 5.0 type 0
part 1 pos 5.9 5.0 5.0 type 0
part 2 pos 6.6 5.5 5.1 type 0

# form a bond of type 0 between particles 1 and 0
part 1 bond 0 0

# form a bond of type 0 between particles 2 and 1
part 2 bond 0 1

# form an angle of type 0 between particles 0 and 2
part 1 bond 1 0 2

integrate 0

puts "\n\n==Program information==\n[code_info]"
puts "\n\n==Interactions==\n[inter]"
puts "\n\n==Energy==\n[analyze energy]"
puts "\n\n==Pressure==\n[analyze pressure]"
puts "\n\n==Stress tensor==\n[analyze stress_tensor]"
puts "\n\n==Particles==\n[part]\n"

integrate 1000

puts "\n\n==Program information==\n[code_info]"
puts "\n\n==Interactions==\n[inter]"
puts "\n\n==Energy==\n[analyze energy]"
puts "\n\n==Pressure==\n[analyze pressure]"
puts "\n\n==Stress tensor==\n[analyze stress_tensor]"
puts "\n\n==Particles==\n[part]\n"
