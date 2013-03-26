import espresso
import logging
from math import sqrt

system, integrator = espresso.standard_system.LennardJones(100, (10,10,10), dt=0.00001, temperature = 1.0)

# logging.getLogger("ExtAnalyze").setLevel(logging.INFO)

print "warming up ..."
capForce = espresso.integrator.CapForce(system, capForce=10000.0)
integrator.addExtension(capForce)
integrator.run(50000)
capForce.disconnect()
print "equilibrating ..."
integrator.dt=0.005
integrator.run(50000)

T = espresso.analysis.Temperature(system)
PressureTensor = espresso.analysis.PressureTensor(system)
# interval between measurements
interval = 10
ExtAnalyzeT = espresso.integrator.ExtAnalyze(T, interval=interval)
integrator.addExtension(ExtAnalyzeT)
ExtAnalyzePressureTensor = espresso.integrator.ExtAnalyze(PressureTensor, interval=interval)
integrator.addExtension(ExtAnalyzePressureTensor)

print "starting integration ... measuring temperature every ", interval, " steps"
integrator.run(10000)

average_T  = ExtAnalyzeT.getAverage()
errorbar_T = sqrt(ExtAnalyzeT.getVariance())

print "average T = ", average_T, " (+- ", errorbar_T, ")"
