//#include "python.hpp"
//#include "NPart.hpp"
#include "NeighborFluctuation.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "boost/unordered_map.hpp"
#include "Cell.hpp"
#include "map"

using namespace espresso;

using namespace std;

namespace espresso {
  namespace analysis {
    using namespace espresso::iterator;

    
    // 
    python::list NeighborFluctuation::computeValues() const {

      real radsq = radius*radius;
      
      // For the moment we don't need the pairs itself. We just need to know what is the 
      // number of neighbors each particle has.
      //boost::unordered_multimap <longint, longint> pairs;
      //boost::unordered_multimap <longint, longint>::iterator pairs_it = pairs.begin();
      
      multiset<longint> neighbNum;
      
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
     	//cout << "pair (" << it->first->id() << ", " << it->second->id() << ")" << endl;
    	if (distsq <= radsq) {
          neighbNum.insert(it->first->id());
          neighbNum.insert(it->second->id());
     	   //pairs_it = pairs.insert(pairs_it, std::make_pair(it->first->id(), it->second->id()));
     	   //pairs_it = pairs.insert(pairs_it, std::make_pair(it->second->id(), it->first->id()));
     	  //std::cout << "added pair (" << it->first->id() << ", " << it->second->id() << ")" << std::endl;
    	}
      }
      
      int n_part = 0;
      for(CellListIterator it(cl); !it.isDone(); ++it) {
    	  //int n = pairs.count(it->id());
    	  int n = neighbNum.count(it->id());
    	  //std::cout << "n of pid " << it->id() << " = " << n << "  pos:"<< it->position() << std::endl;
    	  if (n > 0) {
    		n_part++;
    	    nsum   += n;
    	    nsqsum += n*n;
    	  }
      }
      
      real nave = (1.0*(real)nsum/(real)n_part);
      real naveTot = 0.0;
      
      nave_sq = nave * nave;
      nsq_ave = (1.0*(real)nsqsum/(real)n_part);

      real myE = nsq_ave - nave_sq;
      real E;

      //std::cout << "nave_sq="<<nave_sq<<" nsq_ave="<<nsq_ave<< " myE="<< myE << std::endl;
      //cout << "nave_sq="<<nave_sq<<" nsq_ave="<<nsq_ave<< " myE="<< myE << endl;

      mpi::all_reduce(*system.comm, myE, E, std::plus<real>());
      
      mpi::all_reduce(*system.comm, nave, naveTot, std::plus<real>());
      
      python::list pyList;
      
      pyList.append( E / (real)( system.comm->size() ) );
      
      pyList.append( naveTot / (real)( system.comm->size() ) ) ;

      return pyList;
    }
    
    real NeighborFluctuation::compute() const {
      return -1.0;
    }

    void NeighborFluctuation::registerPython() {
      using namespace espresso::python;
      class_<NeighborFluctuation, bases< Observable > >
        ("analysis_NeighborFluctuation", init< shared_ptr< System > , real >())
        .def("compute", &NeighborFluctuation::computeValues) 
      ;
    }
  }
}
