#include "python.hpp"
#include "NPart.hpp"
#include "NeighborFluctuation.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "boost/unordered_map.hpp"
#include "Cell.hpp"
#include "map"

using namespace espresso;

namespace espresso {
  namespace analysis {
    using namespace espresso::iterator;

    real NeighborFluctuation::compute() const {

      real radsq = radius*radius;
      boost::unordered_multimap <longint, longint> pairs;
      boost::unordered_multimap <longint, longint>::iterator pairs_it = pairs.begin();
      int nsum    = 0;
      int nsqsum  = 0;
      real nsq_ave = 0.0;
      real nave_sq = 0.0;

      System& system = getSystemRef();

      CellList cl = getSystem()->storage->getRealCells();

      for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
    	Real3D d = it->first->position() - it->second->position();
    	real distsq = d.sqr();
    	// std::cout << "distsq: " << distsq << " radsq: " << radsq << std::endl;
    	if (distsq <= radsq) {
     	   pairs_it = pairs.insert(pairs_it, std::make_pair(it->first->id(), it->second->id()));
     	   // std::cout << "added pair (" << it->first->id() << ", " << it->second->id() << ")" << std::endl;
    	}
      }

      int n_part = 0;
      for(CellListIterator it(cl); !it.isDone(); ++it) {
    	  int n = pairs.count(it->id());
    	  // std::cout << "n of pid " << it->id() << " = " << n << std::endl;
    	  if (n > 0) {
    		n_part++;
    	    nsum   += n;
    	    nsqsum += n*n;
    	  }
      }
      nave_sq = (1.0*nsum/n_part) * (1.0*nsum/n_part);
      nsq_ave = (1.0*nsqsum/n_part);

      real myE = nsq_ave - nave_sq;
      real E;

      // std::cout << "nave_sq="<<nave_sq<<" nsq_ave="<<nsq_ave<< " myE="<< myE << std::endl;

      mpi::all_reduce(*system.comm, myE, E, std::plus<real>());

      return E / ( system.comm->size() );
    }

    void NeighborFluctuation::registerPython() {
      using namespace espresso::python;
      class_<NeighborFluctuation, bases< Observable > >
        ("analysis_NeighborFluctuation", init< shared_ptr< System > , real >())
      ;
    }
  }
}
