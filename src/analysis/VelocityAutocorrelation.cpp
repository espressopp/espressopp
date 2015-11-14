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

#include "VelocityAutocorrelation.hpp"

using namespace std;
//using namespace espressopp;

namespace espressopp {
  namespace analysis {

    using namespace iterator;
    
    // calc <vx(0) * vx(t)>
    python::list VelocityAutocorrelation::compute() const{
      
      int M = getListSize();
      real* totZ;
      totZ = new real[M];
      real* Z;
      Z = new real[M];

      python::list pyli;
      
      System& system = getSystemRef();
      
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
      for(int m=0; m<M; m++){
        totZ[m] = 0.0;
        Z[m] = 0.0;
        for(int n=0; n<M-m; n++){
          for (vector<longint>::iterator itr=localIDs.begin(); itr!=localIDs.end(); ++itr) {
            size_t i = *itr;
            Real3D vel1 = getConf(n + m)->getCoordinates(i);
            Real3D vel2 = getConf(n)->getCoordinates(i);
            Z[m] += vel1 * vel2;
          }
        }
        /*
         * additional calculations slow down routine but from the other hand
         * it helps to monitor progress
         */
        if(print_progress && system.comm->rank()==0){
          perc = (int)(m*denom);
          if(perc%5==0){
            cout<<"calculation progress (velocity autocorrelation): "<< perc << " %\r"<<flush;
          }
        }
      }
      if(system.comm->rank()==0)
        cout<<"calculation progress (velocity autocorrelation): 100 %" <<endl;
      
      boost::mpi::all_reduce( *system.comm, Z, M, totZ, plus<real>() );
      
      for(int m=0; m<M; m++){
        totZ[m] /= (real)(M - m);
      }
      
      real inv_coef = 1.0 / (3.0 * num_of_part);
      
      for(int m=0; m<M; m++){
        totZ[m] *= inv_coef;
        pyli.append( totZ[m] );
      }
      
      delete [] Z;
      Z = NULL;
      delete [] totZ;
      totZ = NULL;
      
      return pyli;
    }
    
    // Python wrapping
    void VelocityAutocorrelation::registerPython() {
      using namespace espressopp::python;

      class_<VelocityAutocorrelation, bases<ConfigsParticleDecomp> >
      ("analysis_VelocityAutocorrelation",init< shared_ptr< System > >() )
        .add_property("print_progress", &VelocityAutocorrelation::getPrint_progress, &VelocityAutocorrelation::setPrint_progress)
      ;
    }
  }
}
