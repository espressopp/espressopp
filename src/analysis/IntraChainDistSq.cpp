#include "IntraChainDistSq.hpp"
#include "storage/Storage.hpp"
#include "FixedPairList.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "mpi.h"
#include <map>

using namespace espresso;

namespace espresso {
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
      using namespace espresso::python;

      class_<IntraChainDistSq, bases< AllParticlePos > >
        ("analysis_IntraChainDistSq", init< shared_ptr< System >, shared_ptr< FixedPairList > >())
        .def("compute", &IntraChainDistSq::compute)
        ;

    }


  }
}

