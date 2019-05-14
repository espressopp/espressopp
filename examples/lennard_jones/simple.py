import espressopp
system             = espressopp.System()
integrator     = espressopp.integrator.VelocityVerlet(system)
thermostat             = espressopp.integrator.LangevinThermostat(system)
integrator.addExtension(thermostat)
