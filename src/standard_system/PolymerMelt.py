#  Copyright (C) 2012,2013,2017(H)
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
espressopp.standard_system.PolymerMelt
**************************************


.. function:: espressopp.standard_system.PolymerMelt(num_chains, monomers_per_chain, box, bondlen, rc, skin, dt, epsilon, sigma, shift, temperature, xyzfilename, xyzrfilename)

		:param num_chains: 
		:param monomers_per_chain: 
		:param box: (default: (000))
		:param bondlen: (default: 0.97)
		:param rc: (default: 1.12246)
		:param skin: (default: 0.3)
		:param dt: (default: 0.005)
		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param shift: (default: 'auto')
		:param temperature: (default: None)
		:param xyzfilename: (default: None)
		:param xyzrfilename: (default: None)
		:type num_chains: 
		:type monomers_per_chain: 
		:type box: 
		:type bondlen: real
		:type rc: real
		:type skin: real
		:type dt: real
		:type epsilon: real
		:type sigma: real
		:type shift: 
		:type temperature: 
		:type xyzfilename: 
		:type xyzrfilename: 
		
		returns random walk polymer melt system and integrator:
		if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
"""
import espressopp
import mpi4py.MPI as MPI
import sys

def PolymerMelt(num_chains, monomers_per_chain, box=(0,0,0), bondlen=0.97, rc=1.12246, skin=0.3, dt=0.005, epsilon=1.0, sigma=1.0, shift='auto', temperature=None, xyzfilename=None, xyzrfilename=None):

  if xyzfilename and xyzrfilename:
     print "ERROR: only one of xyzfilename (only xyz data) or xyzrfilename (additional particle radius data) can be provided."
     sys.exit(1)

  if xyzrfilename: 
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf, radiusf = espressopp.tools.readxyzr(xyzrfilename)
    box = (Lxf, Lyf, Lzf)
  elif xyzfilename:
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = espressopp.tools.readxyz(xyzfilename)
    box = (Lxf, Lyf, Lzf)
  else:
    if box[0]<=0 or box[1]<=0 or box[2]<=0:
      print "WARNING: no valid box size specified, box size set to (100,100,100) !"
      box = (100,100,100)

  system         = espressopp.System()
  system.rng     = espressopp.esutil.RNG()
  system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rc,skin)
  cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
  interaction    = espressopp.interaction.VerletListLennardJones(espressopp.VerletList(system, cutoff=rc))
  interaction.setPotential(type1=0, type2=0, potential=espressopp.interaction.LennardJones(epsilon, sigma, rc, shift))
  system.addInteraction(interaction)
  
  integrator     = espressopp.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espressopp.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)

  mass     = 1.0  

  if xyzrfilename: 
    props    = ['id', 'type', 'mass', 'pos', 'v', 'radius']
    bondlist = espressopp.FixedPairList(system.storage)
    for i in xrange(num_chains):
      chain = []
      bonds = []
      for k in xrange(monomers_per_chain):
        idx  =  i * monomers_per_chain + k
        part = [ pidf[idx], typef[idx], mass,
                 espressopp.Real3D(xposf[idx],yposf[idx],zposf[idx]),
                 espressopp.Real3D(xvelf[idx],yvelf[idx],zvelf[idx]),
                 radiusf[idx] ]
        chain.append(part)
        if k>0:
          bonds.append((pidf[idx-1], pidf[idx]))
      system.storage.addParticles(chain, *props)
      system.storage.decompose()
      bondlist.addBonds(bonds)
  elif xyzfilename: 
    props    = ['id', 'type', 'mass', 'pos', 'v']
    bondlist = espressopp.FixedPairList(system.storage)
    for i in xrange(num_chains):
      chain = []
      bonds = []
      for k in xrange(monomers_per_chain):
        idx  =  i * monomers_per_chain + k
        part = [ pidf[idx], typef[idx], mass,
                 espressopp.Real3D(xposf[idx],yposf[idx],zposf[idx]),
                 espressopp.Real3D(xvelf[idx],yvelf[idx],zvelf[idx])]
        chain.append(part)
        if k>0:
          bonds.append((pidf[idx-1], pidf[idx]))
      system.storage.addParticles(chain, *props)
      system.storage.decompose()
      bondlist.addBonds(bonds)
  else:            
    props    = ['id', 'type', 'mass', 'pos', 'v']
    vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)
    bondlist = espressopp.FixedPairList(system.storage)
    pid      = 1
    type     = 0
    chain    = []
    for i in xrange(num_chains):
      startpos = system.bc.getRandomPos()
      positions, bonds = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
      for k in xrange(monomers_per_chain):  
        part = [pid + k, type, mass, positions[k], vel_zero]
        chain.append(part)
      pid += monomers_per_chain
      type += 1
      system.storage.addParticles(chain, *props)
      system.storage.decompose()
      chain = []
      bondlist.addBonds(bonds)

  system.storage.decompose()

  # FENE bonds
  potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
  interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
  system.addInteraction(interFENE)
    
  return system, integrator
  
  
  class KGMelt:
    def __init__(self, num_chains, chain_len):
      self._num_chains    = num_chains
      self._chain_len     = chain_len
      self._num_particles = num_chains * chain_len
      self._density       = 0.8449
      self._L             = pow(self._num_particles / self._density, 1.0/3.0)
      self._box           = (L, L, L)
      self._system        = espress.System()
