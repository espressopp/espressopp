import espresso
import logging
from math import sqrt

system, integrator = espresso.standard_system.LennardJones(1000, (20,20,20), dt=0.00001, temperature = 1.0)

# logging.getLogger("ExtAnalyze").setLevel(logging.INFO)

print "warming up ..."
capForce = espresso.integrator.CapForce(system, capForce=10000.0)
integrator.addExtension(capForce)
integrator.run(50000)
capForce.disconnect()
print "equilibrating ..."
integrator.dt=0.005
integrator.run(50000)

PressureTensor = espresso.analysis.PressureTensor(system)
# interval between measurements
interval = 10
ExtAnalyzePressureTensor = espresso.integrator.ExtAnalyze(PressureTensor, interval=interval)
integrator.addExtension(ExtAnalyzePressureTensor)

print "starting integration ... measuring pressure tensor every ", interval, " steps"
PressureTensor.reset()
integrator.run(10000)

average_PressureTensor = PressureTensor.getAverageValue()

print "average Pressure Tensor = ", average_PressureTensor[:6]
print "          std deviation = ", average_PressureTensor[6:]
print "number of measurements  = ", PressureTensor.getNumberOfMeasurements()
