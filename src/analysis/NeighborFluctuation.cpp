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

#include "NeighborFluctuation.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "boost/unordered_map.hpp"
#include "Cell.hpp"
#include "map"

using namespace espressopp;

using namespace std;

namespace espressopp {
  namespace analysis {
    using namespace espressopp::iterator;

    
    // 
    python::list NeighborFluctuation::computeValues() const {

      real radsq = radius*radius;
      
      // For the moment we don't need the pairs itself. We just need to know what is the 
      // number of neighbors each particle has.
      //boost::unordered_multimap <longint, longint> pairs;
      //boost::unordered_multimap <longint, longint>::iterator pairs_it = pairs.begin();
      
      multiset<longint> neighbNum;
      
      System& system = getSystemRef();

      CellList cl = getSystem()->storage->getRealCells();
      
      longint myN, systemN;
      myN = system.storage->getNRealParticles();
      mpi::all_reduce(*getSystem()->comm, myN, systemN, std::plus<int>());
      
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
      
      
      real nsum    = 0;
      real nsqsum  = 0;
      real nsq_ave = 0.0;
      real nave_sq = 0.0;
      
      real inv_num_part = 1. / (real)systemN;

      for(CellListIterator it(cl); !it.isDone(); ++it) {
    	  //int n = pairs.count(it->id());
    	  int n = neighbNum.count(it->id());
    	  //std::cout << "n of pid " << it->id() << " = " << n << "  pos:"<< it->position() << std::endl;
    	  if (n > 0) {
    	    nsum   += n   * inv_num_part;
    	    nsqsum += n*n * inv_num_part;
    	  }
      }
      
      //cout << "11111nave="<<nsum<<" nsqsum="<<nsqsum<< endl;
      
      //real nave = (1.0*(real)nsum/(real)n_part);
      
      nave_sq = nsum * nsum;
      nsq_ave = nsqsum;

      real myE = nsq_ave - nave_sq;
      real E;

      //std::cout << "nave_sq="<<nave_sq<<" nsq_ave="<<nsq_ave<< " myE="<< myE << std::endl;
      //cout << "nave_sq="<<nave_sq<<" nsq_ave="<<nsq_ave<< " myE="<< myE << endl;

      mpi::all_reduce(*system.comm, myE, E, std::plus<real>());
      
      real naveTot = 0.0;
      mpi::all_reduce(*system.comm, nsum, naveTot, std::plus<real>());
      
      python::list pyList;
      
      pyList.append( E / (real)( system.comm->size() ) );
      
      pyList.append( naveTot / (real)( system.comm->size() ) ) ;

      return pyList;
    }
    
    real NeighborFluctuation::compute() const {
      return -1.0;
    }

    void NeighborFluctuation::registerPython() {
      using namespace espressopp::python;
      class_<NeighborFluctuation, bases< Observable > >
        ("analysis_NeighborFluctuation", init< shared_ptr< System > , real >())
        .def("compute", &NeighborFluctuation::computeValues) 
      ;
    }
  }
}
