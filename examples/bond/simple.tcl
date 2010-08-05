setmd box_l 10.0 10.0 10.0
setmd skin 0.5
setmd time_step 0.01
thermostat off

inter 0 0 lennard-jones 1.0 1.0 2.5 0.0 0.0
inter 0   fene 30.0 1.5

part 0 pos 5.0 5.0 5.0 type 0
part 1 pos 5.9 5.0 5.0 type 0
part 2 pos 6.6 5.5 5.1 type 0

part 3 pos 5.0 5.0 6.0 type 0
part 4 pos 5.9 5.0 6.0 type 0
part 5 pos 6.6 5.5 6.1 type 0

part 1 bond 0 0
part 2 bond 0 1

part 4 bond 0 3
part 5 bond 0 4

integrate 0

puts "\n\nInteractions==\n[inter]"
puts "\n\nEnergy==\n[analyze energy]"
puts "\n\nPressure==\n[analyze pressure]"
#puts "\n\nStress Tensor==\n[analyze stress_tensor]"
#puts "\n\nParticles==\n[part]"

puts "\n\nIntegrate one step\n"

integrate 100

puts "\n\nInteractions==\n[inter]"
puts "\n\nEnergy==\n[analyze energy]"
puts "\n\nPressure==\n[analyze pressure]"
#puts "\n\nStress Tensor==\n[analyze stress_tensor]"
#puts "\n\nParticles==\n[part]"
