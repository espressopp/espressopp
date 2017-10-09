#!/usr/bin/env python2 
#  Copyright (C) 2015-2017(H)
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#####################################################################################################
#                                                                                                   #
#  ESPResSo++ Python script for a Lennard-Jones standard system including pressure tensor analysis  #
#                                                                                                   #
#####################################################################################################

import espressopp
import logging
from math import sqrt

system, integrator = espressopp.standard_system.LennardJones(1000, (20,20,20), dt=0.00001, temperature = 1.0)

# logging.getLogger("ExtAnalyze").setLevel(logging.INFO)

print "warming up ..."
capForce = espressopp.integrator.CapForce(system, capForce=10000.0)
integrator.addExtension(capForce)
integrator.run(50000)
capForce.disconnect()
print "equilibrating ..."
integrator.dt=0.005
integrator.run(50000)

PressureTensor = espressopp.analysis.PressureTensor(system)
# interval between measurements
interval = 10
ExtAnalyzePressureTensor = espressopp.integrator.ExtAnalyze(PressureTensor, interval=interval)
integrator.addExtension(ExtAnalyzePressureTensor)

print "starting integration ... measuring pressure tensor every ", interval, " steps"
PressureTensor.reset()
integrator.run(10000)

average_PressureTensor = PressureTensor.getAverageValue()

print "average Pressure Tensor = ", average_PressureTensor[:6]
print "          std deviation = ", average_PressureTensor[6:]
print "number of measurements  = ", PressureTensor.getNumberOfMeasurements()
