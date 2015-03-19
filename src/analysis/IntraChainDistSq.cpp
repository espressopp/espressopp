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

#include "IntraChainDistSq.hpp"
#include "storage/Storage.hpp"
#include "FixedPairList.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "mpi.h"
#include <map>

using namespace espressopp;

namespace espressopp {
  namespace analysis {

    using namespace iterator;

    LOG4ESPP_LOGGER(IntraChainDistSq::logger, "IntraChainDistSq");

    python::list IntraChainDistSq::compute() {

      python::list R2N;
  	  
      gatherAllPositions();
      
      //for (ConfMap::iterator itr=AllPositions.begin(); itr != AllPositions.end(); ++itr) {
      // 	size_t id = itr->first;
      // 	sBuf p = itr->second;
      // 	R2N.append(python::make_tuple(id, p.r[0], p.r[1], p.r[2]));
      //}

  	  for (FixedPairList::GlobalPairs::const_iterator it=fpl->getGlobalPairs()->begin(); it != fpl->getGlobalPairs()->end(); it++) {
          R2N.append(python::make_tuple(it->first, it->second));
      }

  	  return R2N;
    }

    void IntraChainDistSq::registerPython() {
      using namespace espressopp::python;

      class_<IntraChainDistSq, bases< AllParticlePos > >
        ("analysis_IntraChainDistSq", init< shared_ptr< System >, shared_ptr< FixedPairList > >())
        .def("compute", &IntraChainDistSq::compute)
        ;

    }


  }
}

