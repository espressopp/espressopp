#include "python.hpp"
#include "AllParticlePos.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "mpi.h"
#include <map>
#include <boost/serialization/map.hpp>

using namespace espresso;

namespace espresso {
  namespace analysis {

    using namespace iterator;

    LOG4ESPP_LOGGER(AllParticlePos::logger, "AllParticlePos");

    void AllParticlePos::gatherAllPositions() {

      System& system = getSystemRef();

      int nprocs = system.comm->size();
      int myrank = system.comm->rank();

      for (int rank_i=0; rank_i<nprocs; rank_i++) {

        ConfMap conf;

        if (rank_i == myrank) {
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            int id = cit->id();
            Real3D& pos = cit->position();
            Int3D& img = cit->image();
            Real3D L = system.bc->getBoxL();
            sBuf p;
            for (int i = 0; i < 3; ++i) p.r[i] = pos[i] + img[i] * L[i];
            conf[id] = p;
          }
    	}

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        for (ConfMap::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
        	size_t id = itr->first;
        	sBuf p = itr->second;
            AllPositions[id] = p;
        }
      }

      numParticles = AllPositions.size();
    }

    void AllParticlePos::registerPython() {
      using namespace espresso::python;

      class_<AllParticlePos, boost::noncopyable>("analysis_AllParticlePos", no_init)
        ;
    }


  }
}

