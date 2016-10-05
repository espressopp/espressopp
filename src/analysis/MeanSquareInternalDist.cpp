/*
  Copyright (C) 2012,2013,2016
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

#include "MeanSquareInternalDist.hpp"
#include <boost/serialization/map.hpp>

using namespace std;
//using namespace espressopp;

namespace espressopp {
  namespace analysis {

    //using namespace iterator;
    
    /* 
     * calculates the mean square internal distances
     * 
     * example of how to call it from python:
     * 
     * msid = espressopp.analysis.MeanSquareInternalDist(system,chainlength)
     * msid.gather()            # one can use gatherFromFile instead
     * result = msid.compute()  # computes the msid of the last snapshot added
     *
     * file = open("msid.dat","w")
     *   for i in xrange(chainlength-1):
     *   line = "%d %f %f\n" % (i+1,result[i],result[i]/(i+1)) # n, <r^2>, <r^2>/n
     *   file.write(line)
     * file.close()
     * 
    */

      python::list MeanSquareInternalDist::compute() const{
                         
      python::list pyli;
      
      System& system = getSystemRef();
            
      //creating vector which stores particleIDs for each CPU
      vector<longint> localIDs; // localIDs is a vector with PID of particles in each CPU and it contains whole chains
      for (map<size_t,int>::const_iterator itr=idToCpu.begin(); itr!=idToCpu.end(); ++itr) {
        size_t i = itr->first;
        int whichCPU = itr->second;
        if(system.comm->rank()==whichCPU){
          localIDs.push_back(i); 
        }
      }
      sort(localIDs.begin(), localIDs.end()); //sorts entries from low to high - should not be necessary as long as the above iterator iterates in ascending order according to the keys
                
      // MSID calculation
      
      int M = getListSize(); //number of snapshots/configurations
                 
      real intdist_sum = 0;
                  
      int N_chains = localIDs.size() / chainlength; // N_chains is number of chains in this cpu and chainlength is a global variable from ConfigsParticleDecomp.hpp

      for (int n=1; n<chainlength;n++){ // loop over chainlength

          real intdist = 0.0;
      
          for (int j=0; j<N_chains;j++){
              
              real sumdistsq = 0.0 ;
                                    
              for (int i=j*chainlength; i<(j*chainlength)+chainlength-n;i++){
                  Real3D pos1 = getConf(M-1)->getCoordinates(localIDs[i]);   // compute the MSID of the last added snapshot
                  Real3D pos2 = getConf(M-1)->getCoordinates(localIDs[i+n]); // compute the MSID of the last added snapshot                                                     
                  Real3D delta = pos2 - pos1;   //vector with the distances from i to chainlength-n
                  sumdistsq += delta.sqr();// dx*dx+dy*dy+dz*dz
              }
              
              intdist += sumdistsq/real(chainlength-n);  //sum of the internal distances (chainlength-n) of the chains in this cpu
          }
                    
          boost::mpi::reduce(*system.comm, intdist, intdist_sum, std::plus<real>(), 0);
          
          if (system.comm->rank() == 0) { 
              int chains_total = num_of_part / chainlength;
              pyli.append( intdist_sum / real(chains_total)); 
          }          
      } 
      return pyli;
    }
          
    // Python wrapping
    void MeanSquareInternalDist::registerPython() {
      using namespace espressopp::python;

      class_<MeanSquareInternalDist, bases<ConfigsParticleDecomp> >
      ("analysis_MeanSquareInternalDist", init< shared_ptr< System >, int >() )
      .add_property("print_progress", &MeanSquareInternalDist::getPrint_progress, &MeanSquareInternalDist::setPrint_progress)
      ;
    }
  }
}
