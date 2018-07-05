#  Copyright (C) 2017,2018
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


r"""
**************************************
espressopp.integrator.PIAdressIntegrator
**************************************

The PIAdressIntegrator implements the integration method for Hamiltonian Adaptive Resolution Path Integral Simulations proposed in J. Chem. Phys 147, 244104 (2017) (PI-AdResS). It can be used to run path integral molecular dynamics as well as ring polymer and centroid molecular dynamics in a quantum-classical adaptive resolution fashion, using different empirical force fields. To facilitate an efficient integration, the integrator uses a 3-layer RESPA (J. Chem. Phys. 97, 1990 (1992)) multiple timestepping scheme (inner level: intraatomic spring forces between the Trotter beads. medium level: interatomic bonded forces. outer level: interatomic non-bonded forces). Importantly, the integrator should only be used in combination with PI-AdResS interactions. Furthermore, the integrator has its own thermostat (Langevin), and the only extensions that should be used with it are the Free Energy Compensation (FreeEnergyCompensation) and the Thermodynamic Force (TDforce).

Example:

>>> integrator = espressopp.integrator.PIAdressIntegrator(system, verletlist, timestep_short, timesteps_outerlevel, timesteps_centrallevel, nTrotter, realkinmass, constkinmass, temperature, gamma, centroidthermostat, CMDparameter, PILE, PILElambda, clmassmultiplier, speedupFreezeRings, KTI)
>>> ...
>>> integrator.run(nsteps)

.. function:: espressopp.integrator.PIAdressIntegrator(system, verletlist, timestep, sSteps, mSteps, nTrotter, realKinMass, constKinMass, temperature, gamma, centroidThermostat, CMDparameter, PILE, PILElambda, CLmassmultiplier, speedup, KTI)

        Constructs the PIAdressIntegrator object. Note that all parameters can also be set and fetched via setter and getter functions. Additionally, all parameters except the system and the Verletlist are implemented as class variables that can be directly accessed and modified.

        :param system: system object
        :param verletlist: Verletlist object. Should be an AdResS Verletlist
        :param timestep: (default: 0.0) the inner (shortest) timestep for the calculation of the intraatomic spring forces between the Trotter beads
        :param sSteps: (default: 1) multiplier to construct medium timestep (interatomic bonded forces) as mediumstep = sSteps * timestep
        :param mSteps: (default: 1) multiplier to construct longest timestep (interatomic non-bonded forces) as longstep = mSteps * sSteps * timestep
        :param nTrotter: (default: 32) Trotter number. Should be even and greather than zero.
        :param realKinMass: (default: True) Flag to choose whether to use real kinetic masses. If False, the higher modes' kinetic masses are multiplied with their corresponding eigenvalues of the normal mode transformation. In this way, all higher modes oscillate with the same frequency. If True, we use the kinetic masses for the higher modes which corresponding to the real dynamics (see J. Chem. Phys 147, 244104 (2017) for details)
        :param constKinMass: (default: False) If False, the higher modes' kinetic masses also adaptively change (AKM scheme in J. Chem. Phys 147, 244104 (2017)). If True, the higher modes' kinetic masses are constant throughout the system (CKM scheme in J. Chem. Phys 147, 244104 (2017))
        :param temperature: (default: 2.494353 - this corresponds to 300 Kelvin) the temperature in gromacs units (Boltzmann constant kb is 1)
        :param gamma: (default: 1.0) the Langevin thermostat's friction parameter in 1/ps
        :param centroidThermostat: (default: True) If True, the centroid mode is also thermostated, otherwise only the higher modes' (relevant for centroid molecular dynamics)
        :param CMDparameter: (default: 1.0) The gamma^2 parameter used in centroid molecular dynamics. The higher modes' kinetic masses are rescaled by CMDparameter
        :param PILE: (default: True) If True, the higher modes are thermostated according to the PILE scheme by Ceriotti et al. (J. Chem. Phys 133, 124104 (2010)). Only makes sense in combination when using real kinetic masses (realKinMass = True)
        :param PILElambda: (default: 0.5) lambda parameter to rescale the friction matrix. Default should be good for most applications (J. Chem. Phys 140, 234116 (2014))
        :param CLmassmultiplier: (default: 100.0) multiplier by which the higher modes' spring masses (if constKinMass = False also the kinetic masses) are increased in the classical region
        :param speedup: (default: True) If True, the higher modes' are not integrated in the classical region and also the intraatomistic forces between the Trotter beads are not calculated in the classical region
        :param KTI: (default: False) If True, the particles' resolution parameters and adaptive masses are not updated but can be set by hand everywhere. This is necessary when running Kirkwood Thermodynamic Integration (KTI)
        :type system: shared_ptr<System>
        :type verletlist: shared_ptr<VerletListAdress>
        :type timestep: real
        :type sSteps: int
        :type mSteps: int
        :type nTrotter: int
        :type realKinMass: bool
        :type constKinMass: bool
        :type temperature: real
        :type gamma: real
        :type centroidThermostat: bool
        :type CMDparameter: real
        :type PILE: bool
        :type PILElambda: real
        :type CLmassmultiplier: real
        :type speedup: bool
        :type KTI: bool

.. function:: espressopp.integrator.PIAdressIntegrator.setVerletList(verletlist)

        Sets the VerletList.

        :param verletlist: The VerletListAdress object.
        :type verletlist: espressopp.VerletListAdress

.. function:: espressopp.integrator.PIAdressIntegrator.getVerletList()

        Gets the VerletList.

        :return: the Adress VerletList
        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.integrator.PIAdressIntegrator.setTimeStep(timestep)

        Sets the inner (shortest) timestep.

        :param timestep: the inner timestep
        :type timestep: real

.. function:: espressopp.integrator.PIAdressIntegrator.getTimeStep()

        Gets the inner (shortest) timestep.

        :return: the inner timestep
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.setsStep(sSteps)

        Sets the multiplier to construct medium timestep (interatomic bonded forces) as mediumstep = sSteps * timestep.

        :param sSteps: multiplier to construct medium timestep
        :type sSteps: int

.. function:: espressopp.integrator.PIAdressIntegrator.getsStep()

        Gets the multiplier to construct medium timestep (interatomic bonded forces) as mediumstep = sSteps * timestep.

        :return: multiplier to construct medium timestep
        :rtype: int

.. function:: espressopp.integrator.PIAdressIntegrator.setmStep(mSteps)

        Sets the multiplier to construct longest timestep (interatomic non-bonded forces) as longstep = mSteps * sSteps * timestep.

        :param mSteps: multiplier to construct longest timestep
        :type mSteps: int

.. function:: espressopp.integrator.PIAdressIntegrator.getmStep()

        Gets the multiplier to construct longest timestep (interatomic non-bonded forces) as longstep = mSteps * sSteps * timestep.

        :return: multiplier to construct longest timestep
        :rtype: int

.. function:: espressopp.integrator.PIAdressIntegrator.setNtrotter(nTrotter)

        Sets the Trotter number nTrotter. Should be even and greather than zero. Note that when calling this function, also the normal mode transformation matrix and the eigenvalues are recalculated.

        :param ntrotter: the Trotter number
        :type ntrotter: int

.. function:: espressopp.integrator.PIAdressIntegrator.getNtrotter()

        Gets the Trotter number nTrotter.

        :return: the Trotter number
        :rtype: int

.. function:: espressopp.integrator.PIAdressIntegrator.setRealKinMass(realKinMass)

        Sets the real kinetic mass flag.

        :param realKinMass: the real kinetic mass flag
        :type realKinMass: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getRealKinMass()

        Gets the real kinetic mass flag.

        :return: the real kinetic mass flag
        :rtype: bool

.. function:: espressopp.integrator.PIAdressIntegrator.setConstKinMass(constKinMass)

        Sets the constant kinetic mass flag.

        :param constKinMass: the constant kinetic mass flag
        :type constKinMass: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getConstKinMass()

        Gets the constant kinetic mass flag.

        :return: the constant kinetic mass flag
        :rtype: bool

.. function:: espressopp.integrator.PIAdressIntegrator.setTemperature(temperature)

        Sets the temperature (gromacs units with kb = 1).

        :param temperature: the temperature
        :type temperature: real

.. function:: espressopp.integrator.PIAdressIntegrator.getTemperature()

        Gets the temperature (gromacs units with kb = 1).

        :return: the temperature
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.setGamma(gamma)

        Sets the friction constant gamma (in 1/ps).

        :param gamma: the friction constant gamma
        :type gamma: real

.. function:: espressopp.integrator.PIAdressIntegrator.getGamma()

        Gets the friction constant gamma (in 1/ps).

        :return: the friction constant gamma
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.setCentroidThermostat(centroidThermostat)

        Sets the centroid thermostat flag.

        :param centroidThermostat: the centroid thermostat flag
        :type centroidThermostat: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getCentroidThermostat()

        Gets the centroid thermostat flag.

        :return: the centroid thermostat flag
        :rtype: bool

.. function:: espressopp.integrator.PIAdressIntegrator.setCMDparameter(CMDparameter)

        Sets the centroid molecular dynamics parameter gamma^2 for scaling the kinetic mass.

        :param CMDparameter: the CMD parameter gamma^2
        :type CMDparameter: real

.. function:: espressopp.integrator.PIAdressIntegrator.getCMDparameter()

        Gets the centroid molecular dynamics parameter gamma^2 for scaling the kinetic mass.

        :return: the CMD parameter gamma^2
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.setPILE(PILE)

        Sets the PILE flag.

        :param PILE: the PILE flag
        :type PILE: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getPILE()

        Gets the PILE flag.

        :return: the PILE flag
        :rtype: bool

.. function:: espressopp.integrator.PIAdressIntegrator.setPILElambda(PILElambda)

        Sets the scaling parameter lambda of the PILE thermostat.

        :param PILElambda: the scaling parameter lambda
        :type PILElambda: real

.. function:: espressopp.integrator.PIAdressIntegrator.getPILElambda()

        Gets the scaling parameter lambda of the PILE thermostat.

        :return: the scaling parameter lambda
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.setClmassmultiplier(CLmassmultiplier)

        Sets the multiplier for the higher modes' spring masses in the classical region.

        :param CLmassmultiplier: the classical spring mass multiplier
        :type CLmassmultiplier: real

.. function:: espressopp.integrator.PIAdressIntegrator.getClmassmultiplier()

        Gets the multiplier for the higher modes' spring masses in the classical region.

        :return: the classical spring mass multiplier
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.setSpeedup(speedup)

        Sets the speedup flag.

        :param speedup: the speedup flag
        :type speedup: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getSpeedup()

        Gets the speedup flag.

        :return: the speedup flag
        :rtype: bool

.. function:: espressopp.integrator.PIAdressIntegrator.setKTI(KTI)

        Sets the KTI flag.

        :param speedup: the KTI flag
        :type speedup: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getKTI()

        Gets the KTI flag.

        :return: the KTI flag
        :rtype: bool

.. function:: espressopp.integrator.PIAdressIntegrator.getVerletlistBuilds()

        Gets the number of Verletlist builds.

        :return: number of Verletlist builds
        :rtype: int

.. function:: espressopp.integrator.PIAdressIntegrator.computeRingEnergy()

        Calculates the total configurational energy of all ring polymers in the system based on the springs between the Trotter beads (calculation done using mode coordinates).

        :return: total configurational ring polymer energy
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.computeRingEnergyRaw()

        Calculates the total configurational energy of all ring polymers in the system based on the springs between the Trotter beads (calculation done using the Trotter beads' real space positions).

        :return: total configurational ring polymer energy
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.computeKineticEnergy()

        Calculates the total kinetic energy using the modes' momenta.

        :return: total kinetic energy
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.computePositionDrift(parttype)

        Calculates the average drift force due to the position-dependent spring masses (see Section 5.C. Eq. 63 in J. Chem. Phys 147, 244104 (2017)) on particles of type parttype. To be used during KTI for construction of free energy compensation.

        :param parttype: the particle or atom type
        :type parttype: int
        :return: average drift force due to the position-dependent spring masses
        :rtype: real

.. function:: espressopp.integrator.PIAdressIntegrator.computeMomentumDrift(parttype)

        Calculates the average drift force due to the position-dependent kinetic masses (see Section 5.C. Eq. 62 in J. Chem. Phys 147, 244104 (2017)) on particles of type parttype. To be used during KTI for construction of free energy compensation.

        :param parttype: the particle or atom type
        :type parttype: int
        :return: average drift force due to the position-dependent kinetic masses
        :rtype: real

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.MDIntegrator import *
from _espressopp import integrator_PIAdressIntegrator

import numpy as np

class PIAdressIntegratorLocal(MDIntegratorLocal, integrator_PIAdressIntegrator):
    'The (local) PIAdress Integrator.'
    def __init__(self, system, verletlist, timestep = 0.0, sSteps = 1, mSteps = 1, nTrotter = 32, realKinMass = True, constKinMass = False, temperature = 2.494353, gamma = 1.0, centroidThermostat = True, CMDparameter = 1.0, PILE = True, PILElambda = 0.5, CLmassmultiplier = 100.0, speedup = True, KTI = False):
        if mSteps <= 0 or sSteps <= 0:
            raise ValueError('mSteps and sSteps must be larger than zero. Your inputs: mSteps={}, sSteps={}'.format(mSteps, sSteps))
        if nTrotter <= 0 or nTrotter % 2 != 0:
            raise ValueError('nTrotter must be even and larger than zero. Your input: {}'.format(nTrotter))
        if temperature < 0.0:
            raise ValueError('temperature must be larger or equal zero. Your input: {}'.format(temperature))
        if gamma < 0.0:
            raise ValueError('gamma must be larger or equal zero. Your input: {}'.format(gamma))
        if CMDparameter <= 0.0:
            raise ValueError('CMDparameter must be larger than zero. Your input: {}'.format(CMDparameter))
        if PILElambda < 0.0:
            raise ValueError('PILElambda must be larger or equal zero. Your input: {}'.format(PILElambda))
        if CLmassmultiplier <= 0.0:
            raise ValueError('CLmassmultiplier must be larger than zero. Your input: {}'.format(CLmassmultiplier))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_PIAdressIntegrator, system, verletlist)
            self.setTimeStep(timestep)
            self.setsStep(sSteps)
            self.setmStep(mSteps)
            self.setNtrotter(nTrotter)
            self.setRealKinMass(realKinMass)
            self.setConstKinMass(constKinMass)
            self.setTemperature(temperature)
            self.setGamma(gamma)
            self.setCentroidThermostat(centroidThermostat)
            self.setCMDparameter(CMDparameter)
            self.setPILE(PILE)
            self.setPILElambda(PILElambda)
            self.setClmassmultiplier(CLmassmultiplier)
            self.setSpeedup(speedup)
            self.setKTI(KTI)

    def setTimeStep(self, timestep):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setTimeStep(self, timestep)

    def getTimeStep(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getTimeStep(self)

    def setsStep(self, sSteps):
        if sSteps <= 0:
            raise ValueError('sSteps must be larger than zero. Your input: {}'.format(sSteps))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setsStep(self, sSteps)

    def getsStep(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getsStep(self)

    def setmStep(self, mSteps):
        if mSteps <= 0:
            raise ValueError('mSteps must be larger than zero. Your input: {}'.format(mSteps))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setmStep(self, mSteps)

    def getmStep(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getmStep(self)

    def setNtrotter(self, nTrotter):
        if nTrotter <= 0:
            raise ValueError('nTrotter must be larger than zero. Your input: {}'.format(nTrotter))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setNtrotter(self, nTrotter)

            # Calculate eigenvalues of the normal mode transformation
            for i in range(0,nTrotter):
                self.cxxclass.addEigenValues(self, 4.0*np.sin(np.pi*i/nTrotter)*np.sin(np.pi*i/nTrotter))

            # Calculate the transformation matrix of the normal mode transformation
            transposed_eigenvectors = []
            for k in range(1,nTrotter+1):
                vector = []
                for i in range(0,nTrotter):
                    if i == 0:
                        self.cxxclass.addEVcomponent(self, 1.0/np.sqrt(nTrotter))
                    if i != 0 and i < nTrotter/2:
                        self.cxxclass.addEVcomponent(self, np.sqrt(2.0/nTrotter) * np.cos(2.0*np.pi*k*i/nTrotter))
                    if i == nTrotter/2:
                        self.cxxclass.addEVcomponent(self, (1.0/np.sqrt(nTrotter)) * ((-1.0)**k))
                    if i > nTrotter/2:
                        self.cxxclass.addEVcomponent(self, np.sqrt(2.0/nTrotter) * np.sin(2.0*np.pi*k*i/nTrotter))
                self.cxxclass.addTransposedEigenVector(self)
            self.cxxclass.transp(self)

    def getNtrotter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNtrotter(self)

    def setTemperature(self, temperature):
        if temperature < 0.0:
            raise ValueError('temperature must be larger or equal zero. Your input: {}'.format(temperature))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setTemperature(self, temperature)

    def getTemperature(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getTemperature(self)

    def setGamma(self, gamma):
        if gamma < 0.0:
            raise ValueError('gamma must be larger or equal zero. Your input: {}'.format(gamma))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setGamma(self, gamma)

    def getGamma(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getGamma(self)

    def setCMDparameter(self, CMDparameter):
        if CMDparameter <= 0.0:
            raise ValueError('CMDparameter must be larger than zero. Your input: {}'.format(CMDparameter))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setCMDparameter(self, CMDparameter)

    def getCMDparameter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCMDparameter(self)

    def setPILElambda(self, PILElambda):
        if PILElambda < 0.0:
            raise ValueError('PILElambda must be larger or equal zero. Your input: {}'.format(PILElambda))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPILElambda(self, PILElambda)

    def getPILElambda(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPILElambda(self)

    def setClmassmultiplier(self, CLmassmultiplier):
        if CLmassmultiplier <= 0.0:
            raise ValueError('CLmassmultiplier must be larger than zero. Your input: {}'.format(CLmassmultiplier))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setClmassmultiplier(self, CLmassmultiplier)

    def getClmassmultiplier(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getClmassmultiplier(self)

    def setSpeedup(self, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, speedup)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getSpeedup(self)

    def setKTI(self, KTI):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setKTI(self, KTI)

    def getKTI(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getKTI(self)

    def setCentroidThermostat(self, centroidThermostat):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setCentroidThermostat(self, centroidThermostat)

    def getCentroidThermostat(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCentroidThermostat(self)

    def setPILE(self, PILE):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPILE(self, PILE)

    def getPILE(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPILE(self)

    def setRealKinMass(self, realKinMass):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setRealKinMass(self, realKinMass)

    def getRealKinMass(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getRealKinMass(self)

    def setConstKinMass(self, constKinMass):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setConstKinMass(self, constKinMass)

    def getConstKinMass(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getConstKinMass(self)

    def setVerletList(self, verletlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setVerletList(self, verletlist)

    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

    def computeRingEnergy(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeRingEnergy(self)

    def computeRingEnergyRaw(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeRingEnergyRaw(self)

    def computeKineticEnergy(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeKineticEnergy(self)

    def computePositionDrift(self, parttype):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computePositionDrift(self, parttype)

    def computeMomentumDrift(self, parttype):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeMomentumDrift(self, parttype)

    def getVerletlistBuilds(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletlistBuilds(self)

if pmi.isController :
    class PIAdressIntegrator(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.integrator.PIAdressIntegratorLocal',
          pmiproperty = ['timestep', 'sSteps', 'mSteps', 'nTrotter', 'gamma', 'CMDparameter', 'PILElambda', 'temperature', 'CLmassmultiplier', 'speedup', 'KTI', 'constKinMass', 'verletList', 'centroidThermostat', 'PILE', 'realKinMass', 'verletlistBuilds' ],
          pmicall = ['setTimeStep', 'setmStep', 'setsStep', 'setNtrotter', 'setTemperature', 'setGamma', 'setCMDparameter', 'setPILElambda', 'setClmassmultiplier', 'setSpeedup', 'setKTI', 'setPILE', 'setRealKinMass', 'setCentroidThermostat', 'setConstKinMass', 'setVerletList', 'computeKineticEnergy', 'computeRingEnergy', 'computeRingEnergyRaw', 'computeMomentumDrift', 'computePositionDrift', 'getTimeStep', 'getmStep', 'getsStep', 'getNtrotter', 'getTemperature', 'getGamma', 'getCMDparameter', 'getPILElambda', 'getClmassmultiplier', 'getSpeedup', 'getKTI', 'getPILE', 'getRealKinMass', 'getCentroidThermostat', 'getConstKinMass', 'getVerletList', 'getVerletlistBuilds']
        )
