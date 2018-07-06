#  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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
**************************
espressopp.storage.Storage
**************************

This is the base class for all storage objects.
All derived classes implement at least the following methods:

* `decompose()`

   Send all particles to their corresponding cell/cpu

* `addParticle(pid, pos)`:

   Add a particle to the storage

* `removeParticle(pid)`:

   Remove a particle with id number *pid* from the storage.

   >>> system.storage.removeParticle(4)

   There is an example in *examples* folder

* `getParticle(pid)`:

   Get a particle object.
   This can be used to get specific particle information:

   >>> particle = system.storage.getParticle(15)
   >>> print "Particle ID is       : ", particle.id
   >>> print "Particle position is : ", particle.pos

   you cannot use this particle object to modify particle data.
   You have to use the modifyParticle command for that (see below).

* `addAdrParticle(pid, pos, last_pos)`:

   Add an AdResS Particle to the storage

* `setFixedTuplesAdress(fixed_tuple_list)`:

* `addParticles(particle_list, *properties)`:

   This routine adds particles with certain properties to the storage.

   :param particleList: list of particles (and properties) to be added
   :param properties: property strings

   Each particle in the list must be itself a list where each entry corresponds
   to the property specified in properties.

   Example:

   >>> addParticles([[id, pos, type, ... ], ...], 'id', 'pos', 'type', ...)

* `modifyParticle(pid, property, value, decompose='yes')`

   This routine allows to modify any properties of an already existing particle.

   Example:

   >>> modifyParticle(pid, 'pos', Real3D(new_x, new_y, new_z))

* `removeAllParticles()`:

   This routine removes all particles from the storage.

* 'system':

  The property 'system' returns the System object of the storage.

Examples:

>>> s.storage.addParticles([[1, espressopp.Real3D(3,3,3)], [2, espressopp.Real3D(4,4,4)]],'id','pos')
>>> s.storage.decompose()
>>> s.storage.modifyParticle(15, 'pos', Real3D(new_x, new_y, new_z))



.. function:: espressopp.storage.Storage.addAdrATParticle(pid, \*args)

		:param pid:
		:param \*args:
		:type pid:
		:type \*args:
		:rtype:

.. function:: espressopp.storage.Storage.addParticle(pid, pos)

		:param pid:
		:param pos:
		:type pid:
		:type pos:
		:rtype:

.. function:: espressopp.storage.Storage.addParticles(particleList, \*properties)

		:param particleList:
		:param \*properties:
		:type particleList:
		:type \*properties:
		:rtype:

.. function:: espressopp.storage.Storage.clearSavedPositions()

		:rtype:

.. function:: espressopp.storage.Storage.getParticle(pid)

		:param pid:
		:type pid:
		:rtype:

.. function:: espressopp.storage.Storage.getRealParticleIDs()

		:rtype:

.. function:: espressopp.storage.Storage.modifyParticle(pid, property, value)

		:param pid:
		:param property:
		:param value:
		:type pid:
		:type property:
		:type value:
		:rtype:

.. function:: espressopp.storage.Storage.particleExists(pid)

		:param pid:
		:type pid:
		:rtype:

.. function:: espressopp.storage.Storage.printRealParticles()

		:rtype:

.. function:: espressopp.storage.Storage.removeAllParticles()

		:rtype:

.. function:: espressopp.storage.Storage.removeParticle(pid)

		:param pid:
		:type pid:
		:rtype:

.. function:: espressopp.storage.Storage.restorePositions()

		:rtype:

.. function:: espressopp.storage.Storage.savePositions(idList)

		:param idList:
		:type idList:
		:rtype:

.. function:: espressopp.storage.Storage.setFixedTuplesAdress(fixedtuples)

		:param fixedtuples:
		:type fixedtuples:
"""

from espressopp import pmi
from espressopp import Int3D
import mpi4py.MPI as MPI
import logging
from espressopp import toReal3DFromVector, ParticleLocal, Particle
from espressopp.Exceptions import ParticleDoesNotExistHere

class StorageLocal(object):

    logger = logging.getLogger("Storage")

    def particleExists(self, pid):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        if self.cxxclass.lookupRealParticle(self, pid):
          return True
        else:
          return False

    def addParticle(self, pid, pos):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        self.cxxclass.addParticle(self, pid, toReal3DFromVector(pos))

    def removeParticle(self, pid):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.removeParticle(self, pid)

    def removeAllParticles(self):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.removeAllParticles(self)

    def addAdrATParticle(self, pid, *args):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.addAdrATParticle(
                self, pid, toReal3DFromVector(*args)
                )

    def setFixedTuplesAdress(self, fixedtuples):
        if pmi.workerIsActive():
            self.cxxclass.setFixedTuplesAdress(self, fixedtuples)

    def getParticle(self, pid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return ParticleLocal(pid, self)

    def addParticles(self, particleList, *properties):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():

            index_id          = -1
            index_pos         = -1
            index_modepos     = -1
            index_v           = -1
            index_modemom     = -1
            index_f           = -1
            index_fm          = -1
            index_q           = -1
            index_radius      = -1
            index_fradius     = -1
            index_vradius     = -1
            index_type        = -1
            index_mass        = -1
            index_varmass     = -1
            index_adrAT       = -1 # adress AT particle if 1
            index_lambda_adr  = -1
            index_lambda_adrd = -1
            index_state       = -1
            index_pib         = -1

            last_pos = toReal3DFromVector([-99,-99,-99])

            if properties == None:
                # default properities = (id, pos)
                index_id = 0
                index_pos = 1
                nindex = 2
            else:
                nindex = 0
                for val in properties:
                    if val.lower() == "id": index_id = nindex
                    elif val.lower() == "pos": index_pos = nindex
                    elif val.lower() == "modepos": index_modepos = nindex
                    elif val.lower() == "type": index_type = nindex
                    elif val.lower() == "mass": index_mass = nindex
                    elif val.lower() == "varmass": index_varmass = nindex
                    elif val.lower() == "v": index_v = nindex
                    elif val.lower() == "modemom": index_modemom = nindex
                    elif val.lower() == "f": index_f = nindex
                    elif val.lower() == "fm": index_fm = nindex
                    elif val.lower() == "q": index_q = nindex
                    elif val.lower() == "radius": index_radius = nindex
                    elif val.lower() == "fradius": index_fradius = nindex
                    elif val.lower() == "vradius": index_vradius = nindex
                    elif val.lower() == "adrat": index_adrAT = nindex
                    elif val.lower() == "lambda_adr": index_lambda_adr = nindex
                    elif val.lower() == "lambda_adrd": index_lambda_adrd = nindex
                    elif val.lower() == "state": index_state = nindex
                    elif val.lower() == "pib": index_pib = nindex
                    else: raise SyntaxError("unknown particle property: %s"%val)
                    nindex += 1

            if index_id < 0  : raise "particle property id is mandatory"
            if index_pos < 0 : raise "particle property pos is mandatory"

            # we should check at the begin whether all the particles do not exist.
            doWeAddParticles = True
            for particle in particleList:
              pid = particle[index_id]
              if( self.particleExists(pid) ):
                doWeAddParticles = False
                print "WARNING: Particle ", pid, " already exists"

            if(not doWeAddParticles):
              print 'WARNING: Some particles already exist. The list of particles was not added.'
              return

            for particle in particleList:

                # verify that each particle has enough entries, avoids index errors
                if len(particle) != nindex:
                    raise SyntaxError("particle has %d entries, but %d expected"%(len(particle), nindex))

                id = particle[index_id]
                pos = particle[index_pos]

                if index_adrAT >= 0:
                    if particle[index_adrAT] == 0:
                        #print "%d:  addParticle %d, last_pos=pos %f, %f, %f"%(pmi._MPIcomm.rank,id,pos[0], pos[1], pos[2])
                        storedParticle = self.cxxclass.addParticle(self, id, pos)
                        last_pos = pos
                    else:
                        #print "%d:  addAdrATparticle %d, pos %f, %f, %f, last_pos %f, %f, %f"%(pmi._MPIcomm.rank,id,pos[0],pos[1],pos[2],last_pos[0], last_pos[1], last_pos[2])
                        storedParticle = self.cxxclass.addAdrATParticle(self, id, pos, last_pos)
                else:
                    #print "%d:  addParticle %d, last_pos=pos %f, %f, %f"%(pmi._MPIcomm.rank,id,pos[0], pos[1], pos[2])
                    storedParticle = self.cxxclass.addParticle(self, id, pos)

                if storedParticle != None:
                    self.logger.debug("Processor %d stores particle id = %d"%(pmi.rank, id))
                    self.logger.debug("particle property indexes: id=%i pos=%i type=%i mass=%i v=%i f=%i q=%i radius=%i lambda_adr=%i lambda_adrd=%i state=%i"%(index_id,index_pos,index_type,index_mass,index_v,index_f,index_q,index_radius,index_lambda_adr,index_lambda_adrd,index_state))

                    # only the owner processor writes other properties

                    if index_v >= 0:
                        storedParticle.v = particle[index_v]

                    if index_modepos >= 0:
                        storedParticle.modepos = particle[index_modepos]

                    if index_modemom >= 0:
                        storedParticle.modemom = particle[index_modemom]

                    if index_f >= 0:
                        storedParticle.f = particle[index_f]

                    if index_fm >= 0:
                        storedParticle.fm = particle[index_fm]

                    if index_q >= 0:
                        storedParticle.q = particle[index_q]

                    if index_radius >= 0:
                        storedParticle.radius = particle[index_radius]

                    if index_fradius >= 0:
                        storedParticle.fradius = particle[index_fradius]

                    if index_vradius >= 0:
                        storedParticle.vradius = particle[index_vradius]

                    if index_type >= 0:
                        storedParticle.type = particle[index_type]

                    if index_pib >= 0:
                        storedParticle.pib = particle[index_pib]

                    if index_mass >= 0:
                        storedParticle.mass = particle[index_mass]

                    if index_varmass >= 0:
                        storedParticle.varmass = particle[index_varmass]

                    if index_lambda_adr >= 0:
                        storedParticle.lambda_adr = particle[index_lambda_adr]

                    if index_lambda_adrd >= 0:
                        storedParticle.lambda_adrd = particle[index_lambda_adrd]

                    if index_state >= 0:
                        storedParticle.state = particle[index_state]

    def modifyParticle(self, pid, property, value):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():

          if (self.particleExists(pid)):
              #try:
                #if not particle.isGhost:
                  particle = self.getParticle(pid)
                  self.logger.info("particle pid=%i rank=%i" % (pid, pmi.rank))
                  if   property.lower() == "id"   : raise "particles pid cannot be modified !"
                  elif property.lower() == "pos"  : # alway assume unfolded coordinates
                                                    particle.pos       = value
                                                    particle.imageBox  = Int3D(0, 0, 0)
                  elif property.lower() == "modepos"  : particle.modepos = value
                  elif property.lower() == "img"  : particle.imageBox = value
                  elif property.lower() == "type" : particle.type = value
                  elif property.lower() == "pib" : particle.pib = value
                  elif property.lower() == "mass" : particle.mass = value
                  elif property.lower() == "varmass" : particle.varmass = value
                  elif property.lower() == "v"    : particle.v    = value
                  elif property.lower() == "modemom"  : particle.modemom = value
                  elif property.lower() == "f"    : particle.f    = value
                  elif property.lower() == "q"    : particle.q    = value
                  elif property.lower() == "radius" : particle.radius = value
                  elif property.lower() == "fradius" : particle.fradius = value
                  elif property.lower() == "vradius" : particle.vradius = value
                  elif property.lower() == "lambda_adr" : particle.lambda_adr = value
                  elif property.lower() == "lambda_adrd" : particle.lambda_adrd = value
                  elif property.lower() == "state" : particle.state = value
                  else: raise SyntaxError( 'unknown particle property: %s' % property) # UnknownParticleProperty exception is not implemented
              #except ParticleDoesNotExistHere:
               # self.logger.debug("ParticleDoesNotExistHere pid=% rank=%i" % (pid, pmi.rank))
               # pass
          #else:
           # print "WARNING: Particle ", pid, " does not exist and was not modified"


    def savePositions(self, idList):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        for pid in idList:
          self.cxxclass.savePosition(self, pid)

    def restorePositions(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        self.cxxclass.restorePositions(self)

    def clearSavedPositions(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        self.cxxclass.clearSavedPositions(self)

    def getRealParticleIDs(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.getRealParticleIDs(self)

    def printRealParticles(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        for pid in self.getRealParticleIDs():
          p = ParticleLocal(pid,self)
          print "CPU %-3i ID %-5i TYPE %-3i POS(%8.3f, %8.3f, %8.3f)" % (pmi._MPIcomm.rank, p.id, p.type, p.pos[0], p.pos[1], p.pos[2])

if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "decompose", "addParticles", "setFixedTuplesAdress", "removeAllParticles"],
            pmiproperty = [ "system" ],
            pmiinvoke = ["getRealParticleIDs", "printRealParticles"]
            )

        def particleExists(self, pid):
            return pmi.reduce(pmi.BOR, self.pmiobject, 'particleExists', pid)

        def addParticle(self, pid, pos, checkexist=True):
          if checkexist:
            if( self.particleExists(pid) ):
              print "WARNING: Particle ", pid, " already exists. Therefore it was not added."
              return None
            else:
              pmi.call(self.pmiobject, 'addParticle', pid, pos)
              return Particle(pid, self)
          else:
            pmi.call(self.pmiobject, 'addParticle', pid, pos)

        def removeParticle(self, pid):
            n = pmi.reduce(pmi.SUM, self.pmiobject, 'removeParticle', pid)
            if n == 0:
              print "WARNING: Particle ", pid, " does not exist"
            elif n > 1:
              print "ERROR: Particle ",pid, " did exist more than once !"
              print "       This should never happen !!!"

        def modifyParticle(self, pid, property, value):
          if (self.particleExists(pid)):
            pmi.call(self.pmiobject, 'modifyParticle', pid, property, value)
          else:
            print "WARNING: Particle ", pid, " does not exist and was not modified"


        def addAdrATParticle(self, pid, *args):
          if( self.particleExists(pid) ):
            print "WARNING: Particle ", pid, " already exists. Therefore it was not added."
            return None
          else:
            pmi.call(self.pmiobject, 'addAdrATParticle', pid, *args)
            return Particle(pid, self)

        #def setFixedTuples(self, tuples):
        #    pmi.call(self.pmiobject, 'setFixedTuples', tuples)

        def getParticle(self, pid):
          if( self.particleExists(pid) ):
            return Particle(pid, self)
          else:
            print "WARNING: Particle ", pid, " does not exist"
            return None

        def savePositions(self, idList):
          pmi.call(self.pmiobject, 'clearSavedPositions')
          pmi.call(self.pmiobject, 'savePositions', idList)

        def restorePositions(self):
          pmi.call(self.pmiobject, 'restorePositions')

        def clearSavedPositions(self):
          pmi.call(self.pmiobject, 'clearSavedPositions')

