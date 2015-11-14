/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"
#include "AllParticlePos.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "mpi.h"
#include <map>
#include <boost/serialization/map.hpp>

using namespace espressopp;

namespace espressopp {
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
    
    ConfMap AllParticlePos::getAllPositions(){
      return AllPositions;
    }
    

    void AllParticlePos::registerPython() {
      using namespace espressopp::python;

      class_<AllParticlePos, boost::noncopyable>("analysis_AllParticlePos", no_init)
      
        .def("gatherAllPositions", &AllParticlePos::gatherAllPositions)
        ;
    }


  }
}

