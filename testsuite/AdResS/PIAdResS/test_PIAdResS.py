#!/usr/bin/env python
#
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
#
# -*- coding: utf-8 -*-

import espressopp
import mpi4py.MPI as MPI
from espressopp import Real3D
from espressopp.tools import decomp

import unittest

class ParametrizedTestCase(unittest.TestCase):
    """ TestCase classes that want to be parametrized should inherit from this class. """
    def __init__(self, methodName='runTest', param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_class, param=None):
        """ Create a suite containing all tests taken from the given subclass, passing them the parameter 'param'. """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_class)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_class(name, param=param))
        return suite

class TestPIAdResS(ParametrizedTestCase):
    def test_PIAdResS(self):
        print 'param =', self.param

        constkinmass = self.param['constkinmass']
        PILE = self.param['PILE']
        realkinmass = self.param['realkinmass']
        centroidthermostat = self.param['centroidthermostat']
        KTI = self.param['KTI']
        speedupInterAtom = self.param['speedupInterAtom']
        speedupFreezeRings = self.param['speedupFreezeRings']
        spherical_adress = self.param['spherical_adress']
        nb_forcefield_setup = self.param['nb_forcefield_setup']
        energy_before = self.param['energy_before']
        energy_after = self.param['energy_after']

        steps = 10
        timestep_short = 0.001/16.0
        multiplier_short_to_medium = 4
        multiplier_medium_to_long = 4
        interaction_cutoff = 0.84
        potential_cutoff = 0.78
        skin = 0.1
        gamma = 0.0
        temp = 2.50266751
        if KTI:
            ex_size = 100.0
        else:
            ex_size = 1.0
        hy_size = 1.5
        nTrotter = 4
        clmassmultiplier = 100.0
        PILElambda = 0.0
        CMDparameter = 1.0

        tabFEC_H = "FEC_H.dat"
        tabFEC_O = "FEC_O.dat"
        tabTHDF_H = "ThdForce_H.dat"
        tabTHDF_O = "ThdForce_O.dat"
        tabAngle = "tableESP_angle.dat"
        tabBondHH = "tableESP_bondHH.dat"
        tabBondOH = "tableESP_bondOH.dat"
        tabHW_HW = "tableESP_HW_HW.dat"
        tabHW_OW = "tableESP_HW_OW.dat"
        tabOW_OW = "tableESP_OW_OW.dat"

        pid, types, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("input_test.xyz")
        masses = []
        for item in types:
            if item == 1:
                masses.append(15.9994)
            else:
                masses.append(1.008)

        num_Trotter_beads = len(x)
        num_atoms = len(x)/nTrotter
        size = (Lx, Ly, Lz)

        system = espressopp.System()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
        system.skin = skin
        comm = MPI.COMM_WORLD
        nodeGrid = decomp.nodeGrid(comm.size)
        cellGrid = decomp.cellGrid(size, nodeGrid, interaction_cutoff, skin)
        system.rng = espressopp.esutil.RNG()
        system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

        props = ['id', 'pos', 'v', 'f', 'pib', 'type', 'mass', 'adrat']
        allParticlesAT = []
        allParticles = []
        tuples = []

        for pid_trotter in range(num_Trotter_beads):
            allParticlesAT.append([pid_trotter + 1, Real3D(x[pid_trotter], y[pid_trotter], z[pid_trotter]),
                                Real3D(vx[pid_trotter], vy[pid_trotter], vz[pid_trotter]), Real3D(0, 0, 0),
                                pid_trotter%nTrotter + 1, types[pid_trotter], masses[pid_trotter], 1])
        for pid_atom in range(num_atoms):
            tmptuple = [pid_atom+num_Trotter_beads+1]
            for pid_trotter in range(nTrotter):
                pid = pid_atom*nTrotter+pid_trotter
                tmptuple.append((allParticlesAT[pid])[0])
            firstParticleId=tmptuple[1]
            cmp=allParticlesAT[firstParticleId-1][1]
            cmv=allParticlesAT[firstParticleId-1][2]
            allParticles.append([pid_atom+num_Trotter_beads+1, Real3D(cmp[0], cmp[1], cmp[2]), Real3D(cmv[0], cmv[1], cmv[2]), Real3D(0, 0, 0), 0,
                                types[pid_atom*nTrotter], masses[pid_atom*nTrotter], 0])
            for pid_trotter in range(nTrotter):
                pid = pid_atom*nTrotter+pid_trotter
                allParticles.append([(allParticlesAT[pid])[0],
                                    (allParticlesAT[pid])[1], (allParticlesAT[pid])[2], (allParticlesAT[pid])[3], (allParticlesAT[pid])[4],
                                    (allParticlesAT[pid])[5], (allParticlesAT[pid])[6], (allParticlesAT[pid])[7]])
            tuples.append(tmptuple)

        system.storage.addParticles(allParticles, *props)
        ftpl = espressopp.FixedTupleListAdress(system.storage)
        ftpl.addTuples(tuples)
        system.storage.setFixedTuplesAdress(ftpl)
        system.storage.decompose()

        bondsOH = []
        bondsHH = []
        for part in range(num_atoms/3):
            bondsOH.append((num_Trotter_beads + 1 + 3*part, num_Trotter_beads + 1 + 3*part+1))
            bondsOH.append((num_Trotter_beads + 1 + 3*part, num_Trotter_beads + 1 + 3*part+2))
            bondsHH.append((num_Trotter_beads + 1 + 3*part+1, num_Trotter_beads + 1 + 3*part+2))
        fplOH = espressopp.FixedPairList(system.storage)
        fplHH = espressopp.FixedPairList(system.storage)
        fplOH.addBonds(bondsOH)
        fplHH.addBonds(bondsHH)

        angles = []
        for part in range(num_atoms/3):
            angles.append((num_Trotter_beads + 1 + 3*part+1, num_Trotter_beads + 1 + 3*part, num_Trotter_beads + 1 + 3*part+2))
        ftl = espressopp.FixedTripleList(system.storage)
        ftl.addTriples(angles)

        vl = espressopp.VerletListAdress(system, cutoff=interaction_cutoff, adrcut=interaction_cutoff, dEx=ex_size, dHy=hy_size,
                                         adrCenter=[Lx/2, Ly/2, Lz/2], exclusionlist=bondsOH+bondsHH, sphereAdr=spherical_adress)

        if nb_forcefield_setup == 1:
            interNB = espressopp.interaction.VerletListPIadressTabulatedLJ(vl, ftpl, nTrotter, speedupInterAtom)
        elif nb_forcefield_setup == 2:
            interNB = espressopp.interaction.VerletListPIadressNoDriftTabulated(vl, ftpl, nTrotter, speedupInterAtom)
        elif nb_forcefield_setup == 3:
            interNB = espressopp.interaction.VerletListPIadressTabulated(vl, ftpl, nTrotter, speedupInterAtom)
        else:
            raise ValueError("Wrong nb_forcefield_setup integer (only 1,2,3 accepted.")
        potOOqm = espressopp.interaction.Tabulated(itype=3, filename=tabOW_OW, cutoff=potential_cutoff)
        potHOqm = espressopp.interaction.Tabulated(itype=3, filename=tabHW_OW, cutoff=potential_cutoff)
        potHHqm = espressopp.interaction.Tabulated(itype=3, filename=tabHW_HW, cutoff=potential_cutoff)
        if nb_forcefield_setup == 1:
            interNB.setPotentialQM(type1=1, type2=1, potential=potOOqm)
            interNB.setPotentialQM(type1=1, type2=0, potential=potHOqm)
            interNB.setPotentialQM(type1=0, type2=0, potential=potHHqm)
            potOOcl = espressopp.interaction.LennardJones(epsilon=temp, sigma=0.25, shift='auto', cutoff=1.122462048309373*0.25)
            interNB.setPotentialCL(type1=1, type2=1, potential=potOOcl)
        elif nb_forcefield_setup == 2:
            interNB.setPotential(type1=1, type2=1, potential=potOOqm)
            interNB.setPotential(type1=1, type2=0, potential=potHOqm)
            interNB.setPotential(type1=0, type2=0, potential=potHHqm)
        elif nb_forcefield_setup == 3:
            interNB.setPotentialQM(type1=1, type2=1, potential=potOOqm)
            interNB.setPotentialQM(type1=1, type2=0, potential=potHOqm)
            interNB.setPotentialQM(type1=0, type2=0, potential=potHHqm)
            interNB.setPotentialCL(type1=1, type2=1, potential=potOOqm)
            interNB.setPotentialCL(type1=1, type2=0, potential=potHOqm)
            interNB.setPotentialCL(type1=0, type2=0, potential=potHHqm)
        system.addInteraction(interNB)

        potBondHH = espressopp.interaction.Tabulated(itype=3, filename=tabBondHH)
        potBondOH = espressopp.interaction.Tabulated(itype=3, filename=tabBondOH)
        interBondedHH = espressopp.interaction.FixedPairListPIadressTabulated(system, fplHH, ftpl, potBondHH, nTrotter, speedupInterAtom)
        interBondedOH = espressopp.interaction.FixedPairListPIadressTabulated(system, fplOH, ftpl, potBondOH, nTrotter, speedupInterAtom)
        system.addInteraction(interBondedHH)
        system.addInteraction(interBondedOH)

        potAngle = espressopp.interaction.TabulatedAngular(itype=3, filename=tabAngle)
        interAngle = espressopp.interaction.FixedTripleListPIadressTabulatedAngular(system, ftl, ftpl, potAngle, nTrotter, speedupInterAtom)
        system.addInteraction(interAngle)

        integrator = espressopp.integrator.PIAdressIntegrator(system=system, verletlist=vl, timestep=timestep_short, sSteps=multiplier_short_to_medium, mSteps=multiplier_medium_to_long, nTrotter=nTrotter, realKinMass=realkinmass, constKinMass=constkinmass, temperature=temp, gamma=gamma, centroidThermostat=centroidthermostat, CMDparameter=CMDparameter, PILE=PILE, PILElambda=PILElambda, CLmassmultiplier=clmassmultiplier, speedup=speedupFreezeRings, KTI=KTI)

        if not KTI:
            fec = espressopp.integrator.FreeEnergyCompensation(system, center=[Lx/2, Ly/2, Lz/2], ntrotter=nTrotter)
            fec.addForce(itype=3, filename=tabFEC_O, type=1)
            fec.addForce(itype=3, filename=tabFEC_H, type=0)
            integrator.addExtension(fec)
            thdf = espressopp.integrator.TDforce(system, vl)
            thdf.addForce(itype=3, filename=tabTHDF_O, type=1)
            thdf.addForce(itype=3, filename=tabTHDF_H, type=0)
            integrator.addExtension(thdf)

        if KTI:
            for i in range(1, num_Trotter_beads + num_atoms + 1):
                system.storage.modifyParticle(i, 'lambda_adrd', 0.0)
                system.storage.modifyParticle(i, 'lambda_adr', 0.0)
                system.storage.modifyParticle(i, 'varmass', clmassmultiplier*system.storage.getParticle(i).mass)
            system.storage.decompose()

        espressopp.tools.AdressDecomp(system, integrator)

        Eb = interBondedOH.computeEnergy() + interBondedHH.computeEnergy()
        EAng = interAngle.computeEnergy()
        ELj= interNB.computeEnergy()
        Ek = integrator.computeKineticEnergy()
        EPI = integrator.computeRingEnergy()
        if KTI==False:
            Ecorr = fec.computeCompEnergy() + thdf.computeTDEnergy()
        else:
            Ecorr = 0.0
        energy_before_thistest =  Ek+Eb+EAng+ELj+EPI+Ecorr

        integrator.run(steps)

        Eb = interBondedOH.computeEnergy() + interBondedHH.computeEnergy()
        EAng = interAngle.computeEnergy()
        ELj= interNB.computeEnergy()
        Ek = integrator.computeKineticEnergy()
        EPI = integrator.computeRingEnergy()
        if KTI==False:
            Ecorr = fec.computeCompEnergy() + thdf.computeTDEnergy()
        else:
            Ecorr = 0.0
        energy_after_thistest = Ek+Eb+EAng+ELj+EPI+Ecorr

        self.assertAlmostEqual(energy_before_thistest, energy_before, places=5)
        self.assertAlmostEqual(energy_after_thistest, energy_after, places=5)


if __name__ == '__main__':
        suite = unittest.TestSuite()

        # Make reference parameter dictionary
        parameter_dict_ref = {}
        parameter_dict_ref['constkinmass'] = False
        parameter_dict_ref['PILE'] = True
        parameter_dict_ref['realkinmass'] = True
        parameter_dict_ref['centroidthermostat'] = True
        parameter_dict_ref['KTI'] = False
        parameter_dict_ref['speedupInterAtom'] = True
        parameter_dict_ref['speedupFreezeRings'] = False
        parameter_dict_ref['spherical_adress'] = False
        parameter_dict_ref['nb_forcefield_setup'] = 1
        parameter_dict_ref['energy_before'] = 0.0
        parameter_dict_ref['energy_after'] = 0.0

        # Make test cases
        parameter_dict_test1 = dict.copy(parameter_dict_ref)
        parameter_dict_test1['energy_before'] = 2182.05339049
        parameter_dict_test1['energy_after'] = 2181.96191713
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test1))

        parameter_dict_test2 = dict.copy(parameter_dict_ref)
        parameter_dict_test2['energy_before'] = 2161.92463502
        parameter_dict_test2['energy_after'] = 2161.76667819
        parameter_dict_test2['nb_forcefield_setup'] = 2
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test2))

        parameter_dict_test3 = dict.copy(parameter_dict_ref)
        parameter_dict_test3['energy_before'] = 2161.92463502
        parameter_dict_test3['energy_after'] = 2161.76667819
        parameter_dict_test3['nb_forcefield_setup'] = 3
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test3))

        parameter_dict_test4 = dict.copy(parameter_dict_ref)
        parameter_dict_test4['energy_before'] = 1713.79167547
        parameter_dict_test4['energy_after'] = 1713.12411385
        parameter_dict_test4['constkinmass'] = True
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test4))

        parameter_dict_test5 = dict.copy(parameter_dict_ref)
        parameter_dict_test5['energy_before'] = 2182.05339049
        parameter_dict_test5['energy_after'] = 2181.96191713
        parameter_dict_test5['PILE'] = False
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test5))

        parameter_dict_test6 = dict.copy(parameter_dict_ref)
        parameter_dict_test6['energy_before'] = 3277.87198138
        parameter_dict_test6['energy_after'] = 3277.72611308
        parameter_dict_test6['realkinmass'] = False
        parameter_dict_test6['PILE'] = False
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test6))

        parameter_dict_test7 = dict.copy(parameter_dict_ref)
        parameter_dict_test7['energy_before'] = 2182.05339049
        parameter_dict_test7['energy_after'] = 2181.96191713
        parameter_dict_test7['centroidthermostat'] = False
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test7))

        parameter_dict_test8 = dict.copy(parameter_dict_ref)
        parameter_dict_test8['energy_before'] = 29760.7680052
        parameter_dict_test8['energy_after'] = 29760.8107422
        parameter_dict_test8['KTI'] = True
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test8))

        parameter_dict_test9 = dict.copy(parameter_dict_ref)
        parameter_dict_test9['energy_before'] = 1989.29552038
        parameter_dict_test9['energy_after'] = 1989.19640208
        parameter_dict_test9['speedupFreezeRings'] = True
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test9))

        parameter_dict_test10 = dict.copy(parameter_dict_ref)
        parameter_dict_test10['energy_before'] = 2186.59523062
        parameter_dict_test10['energy_after'] = 2186.46702152
        parameter_dict_test10['speedupInterAtom'] = False
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test10))

        parameter_dict_test11 = dict.copy(parameter_dict_ref)
        parameter_dict_test11['energy_before'] = 3637.02668470
        parameter_dict_test11['energy_after'] = 3636.26448522
        parameter_dict_test11['spherical_adress'] = True
        suite.addTest(ParametrizedTestCase.parametrize(TestPIAdResS, param=parameter_dict_test11))

        # Run tests
        unittest.TextTestRunner().run(suite)
