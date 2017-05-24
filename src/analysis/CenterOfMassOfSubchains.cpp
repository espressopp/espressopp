/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research
  
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

#include "CenterOfMassOfSubchains.hpp"
#include "bc/BC.hpp"

using namespace std;
//using namespace espressopp;

namespace espressopp {
  namespace analysis {

    using namespace iterator;
    
    python::list CenterOfMassOfSubchains::compute() const{
      
      int M = getListSize();
      long unsigned int num_of_subchain = num_of_part/chainlength;
      //should be initialized?
      Real3D* totCOM;
      totCOM = new Real3D[num_of_subchain];
      Real3D* COM;
      COM = new Real3D[num_of_subchain];

      python::list pyli;
      
      System& system = getSystemRef();
      Real3D Li = system.bc->getBoxL();

      for(long unsigned int id = 0; id < num_of_subchain; id++) {
	COM[id] = 0.;
	totCOM[id] = 0.;
      }

      //cout << "In compute, the size of idToCPU = " << idToCpu.size() << endl;
      vector<longint> localIDs;
      for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
        size_t i = itr->first;
        int whichCPU = itr->second;
        if(system.comm->rank()==whichCPU){
	    localIDs.push_back(i);
        }
      }
      
      int perc=0;
      real denom = 100.0 / (real)M;
      Real3D cp;
      for(int m = 0; m < M; m++){
	for (vector<longint>::iterator itr=localIDs.begin(); itr!=localIDs.end(); ++itr) {
	  size_t i = *itr;
          Real3D pos_i = getConf(m)->getCoordinates(i);
	  if ((i - 1) % chainlength != 0) {
	    Real3D diff;
	    system.bc->getMinimumImageVectorBox(diff,
						pos_i,
						cp);
	    pos_i = cp + diff;
	  }
	  cp = pos_i;
	  int cid = idToCid.at(i - 1);
	  COM[cid] += pos_i;
	}
        /*
         * additional calculations slow down routine but from the other hand
         * it helps to monitor progress
         */
        if(print_progress && system.comm->rank()==0){
          perc = (int)(m*denom);
          if(perc%5==0){
            cout<<"calculation progress (center of mass): "<< perc << " %\r"<<flush;
          }
        }
      }
      if(print_progress && system.comm->rank()==0)
	cout<<"calculation progress (center of mass): 100 %" <<endl;

      for(long unsigned int id = 0; id < num_of_subchain; id++) {
	COM[id] /= M*chainlength;
	COM[id] = Real3D(fmod(COM[id][0] + Li[0], Li[0]),
			 fmod(COM[id][1] + Li[1], Li[1]),
			 fmod(COM[id][2] + Li[2], Li[2]));
      }
      
      boost::mpi::all_reduce( *system.comm, COM, num_of_subchain, totCOM, plus<Real3D>() );

      for(long unsigned int id = 0; id < num_of_subchain; id++) {
        pyli.append( totCOM[id] );
      }

      delete [] COM;
      COM = NULL;
      delete [] totCOM;
      totCOM = NULL;

      return pyli;
    }
    
    // Python wrapping
    void CenterOfMassOfSubchains::registerPython() {
      using namespace espressopp::python;

      class_<CenterOfMassOfSubchains, bases<ConfigsParticleDecomp> >
	  ("analysis_CenterOfMassOfSubchains",init< shared_ptr< System >, int>() )
        .add_property("print_progress", &CenterOfMassOfSubchains::getPrint_progress, &CenterOfMassOfSubchains::setPrint_progress)
	.def("getPrint_progress", &CenterOfMassOfSubchains::getPrint_progress)
	.def("setPrint_progress", &CenterOfMassOfSubchains::setPrint_progress)
      ;
    }
  }
}
