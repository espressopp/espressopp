from espresso import pmi
import MPI
import logging
from espresso import toReal3DFromVector, ParticleLocal, Particle

class StorageLocal(object):
    """Abstract local base class for storing particles"""

    logger = logging.getLogger("Storage")

    def addParticle(self, pid, *args):
        """Add a particle locally if it belongs to the local domain."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.addParticle(
                self, pid, toReal3DFromVector(*args)
                )
    
    def getParticle(self, pid):
        """Get the local particle. If it is not on this node, any
        attempt to access the particle will raise an exception."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return ParticleLocal(pid, self)

    def addParticles(self, particleList, *properties):
        """
        This routine adds particles with certain properties to the storage
        where only one processor adds the particle in its local storage.

        :param particleList: list of list particles (and properties) to be added
        :param properties: property strings

        Each particle in the list must be itself a list where each entry corresponds
        to the property specified in properties.
        
        Example: addParticles([[id, pos, type, ... ], ...], 'id', 'pos', 'type', ...)

        Improvement: make list of supported properties more flexible in use
        """
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():

            index_id   = -1
            index_pos  = -1
            index_v    = -1
            index_f    = -1
            index_q    = -1
            index_type = -1
            index_mass = -1

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
                    elif val.lower() == "type": index_type = nindex
                    elif val.lower() == "mass": index_mass = nindex
                    elif val.lower() == "v": index_v = nindex
                    elif val.lower() == "f": index_f = nindex
                    elif val.lower() == "q": index_q = nindex
                    else: raise "unknown particle property: %s"%val
                    nindex += 1

            if index_id < 0  : raise "particle property id is mandatory"
            if index_pos < 0 : raise "particle property pos is mandatory"

            for particle in particleList:

                # verify that each particle has enough entries, avoids index errors
                if len(particle) != nindex:
                    raise "particle has %d entries, but %d expected"%(len(particle), nindex)

                id = particle[index_id]
                pos = particle[index_pos]

                storedParticle = self.cxxclass.addParticle(self, id, pos)

                if storedParticle != None:

                    self.logger.debug("Processor %d stores particle id = %d"%(pmi._MPIcomm.rank, id))
                    self.logger.debug("particle property indexes: id=%i pos=%i type=%i mass=%i v=%i f=%i q=%i"%(index_id,index_pos,index_type,index_mass,index_v,index_f,index_q))

                    # only the owner processor writes other properties

                    if index_v >= 0:
                        storedParticle.v = particle[index_v]

                    if index_f >= 0:
                        storedParticle.f = particle[index_f] 

                    if index_q >= 0:
                        storedParticle.q = particle[index_q]

                    if index_type >= 0:
                        storedParticle.type = particle[index_type]

                    if index_mass >= 0:
                        storedParticle.mass = particle[index_mass]
 
    

if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "decompose", "addParticles" ],
            pmiproperty = [ "system" ]
            )

        def addParticle(self, pid, *args):
            pmi.call(self.pmiobject, 'addParticle', pid, *args)
            return Particle(pid, self)

        def getParticle(self, pid):
            return Particle(pid, self)
