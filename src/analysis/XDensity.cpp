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
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "Configuration.hpp"
#include "XDensity.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"
#include "math.h"

#include <boost/serialization/map.hpp>

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {
    
    // splitN is a level of discretisation of density profile (how many elements it contains)
    python::list XDensity::computeArray(int splitN) const {

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      real Li = system.bc->getBoxL()[0];
      //int nprocs = system.comm->size();
      //int myrank = system.comm->rank();
      real * histogram = 0;
      histogram = new real[splitN];
      for(int i=0;i<splitN;i++) histogram[i]=0.0;
      //std::vector<real> histogram;
      //for(int i=0;i<splitN;i++) histogram.push_back(0.0);
      real dr = Li / (real)splitN;
      int num_part = 0;
      //ConfigurationPtr config = make_shared<Configuration> ();
      /*for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, Real3D > conf;
        if (rank_i == myrank) {
            if(system.storage->getFixedTuples()){
                shared_ptr<FixedTupleListAdress> fixedtupleList=system.storage->getFixedTuples();
                CellList realCells = system.storage->getRealCells();
                                
                for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {  // Iterate over all (CG) particles.              
                    Particle &vp = *cit;
                    FixedTupleListAdress::iterator it2;
                    it2 = fixedtupleList->find(&vp);

                    if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
                          std::vector<Particle*> atList;
                          atList = it2->second;
                          for (std::vector<Particle*>::iterator it3 = atList.begin();
                                               it3 != atList.end(); ++it3) {
                              Particle &at = **it3;
                              int id = at.id();
                              conf[id] = at.position();
                          }  
                    }

                    else{   // If not, use CG particle itself for calculation.
                          int id = cit->id();
                          conf[id] = cit->position();
                    }

                }
            
            }
            else{
                CellList realCells = system.storage->getRealCells();
                for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
                  int id = cit->id();
                  conf[id] = cit->position();
                }
            }
    	}

        boost::mpi::broadcast(*system.comm, conf, rank_i);

        // for simplicity we will number the particles from 0
        for (map<size_t,Real3D>::iterator itr=conf.begin(); itr != conf.end(); ++itr) {
          //size_t id = itr->first;
          Real3D p = itr->second;
          config->set(num_part, p[0], p[1], p[2]);
          num_part ++;
        }
      }*/
      // now all CPUs have all particle coords and num_part is the total number of particles
      
      // use all cpus
      // TODO it could be a problem if   n_nodes > num_part
      /*int numi = num_part / nprocs + 1;
      int mini = myrank * numi;
      int maxi = mini + numi;
      
      if(maxi>num_part) maxi = num_part;

      int perc=0, perc1=0;*/

      real pos = 0.0;
      int bin = 0;
      if(system.storage->getFixedTuples()){
            shared_ptr<FixedTupleListAdress> fixedtupleList=system.storage->getFixedTuples();
            CellList realCells = system.storage->getRealCells();

            for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {  // Iterate over all (CG) particles.              
                Particle &vp = *cit;
                FixedTupleListAdress::iterator it2;
                it2 = fixedtupleList->find(&vp);

                if (it2 != fixedtupleList->end()) {  // Are there atomistic particles for given CG particle? If yes, use those for calculation.
                      std::vector<Particle*> atList;
                      atList = it2->second;
                      for (std::vector<Particle*>::iterator it3 = atList.begin();
                                           it3 != atList.end(); ++it3) {
                          Particle &at = **it3;
                          pos = at.position()[0];
                          if (pos < 0.0)
                          {        
                                    bin = floor ((pos+Li)/dr);
                          }
                          else if(pos > Li)
                          {        
                                    bin = floor ((pos-Li)/dr);
                          }
                          else
                          {
                                    bin = floor (pos/dr);
                          }
                          histogram[bin] += 1.0;
                          num_part += 1;
                      }  
                }

                else{   // If not, use CG particle itself for calculation.
                      pos = cit->position()[0];
                      if (pos < 0.0)
                      {        
                                bin = floor ((pos+Li)/dr);
                      }
                      else if(pos > Li)
                      {        
                                bin = floor ((pos-Li)/dr);
                      }
                      else
                      {
                                bin = floor (pos/dr);
                      }
                      histogram[bin] += 1.0;
                      num_part += 1;
                }

            }

      }
      else{
            CellList realCells = system.storage->getRealCells();
            for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
              pos = cit->position()[0];
              if (pos < 0.0)
              {        
                        bin = floor ((pos+Li)/dr);
              }
              else if(pos > Li)
              {        
                        bin = floor ((pos-Li)/dr);
              }
              else
              {
                        bin = floor (pos/dr);
              }
              histogram[bin] += 1.0;
              num_part += 1;
            }
      }

      //int bin = 0;
      /*for(int i = mini; i<maxi; i++){
        Real3D coordP1 = config->getCoordinates(i);
        //int bin = floor (coordP1[0]/dr);
        if (coordP1[0] < 0.0)
        {        
                bin = floor ((coordP1[0]+Li[0])/dr);
                //cout << "bin: " << bin << " \n";
        }
        else if(coordP1[0] > Li[0])
        {        
                bin = floor ((coordP1[0]-Li[0])/dr);
                //cout << "bin: " << bin << " \n";
        }
        else
        {
                bin = floor (coordP1[0]/dr);
                //cout << "bin: " << bin << " \n";
        }
        histogram[bin] += 1.0;
        */
        /*if(system.comm->rank()==0){
          perc = (int)(100*(real)(i-mini)/(real)(maxi-mini));
          if(perc>perc1){
            cout<<"calculation progress (density profile along x axis): "<< perc << " %\r"<<flush;
            perc1 = perc;
          }
        }*/
      //}
      //if(system.comm->rank()==0)
        //cout<<"calculation progress (density profile along x axis): "<< 100 << " %" <<endl;

      int total_num_part = 0;
      //std::vector<real> totHistogram;
      //for(int i=0;i<splitN;i++) totHistogram.push_back(0.0);
      real * totHistogram = 0;
      totHistogram = new real[splitN];
      for(int i=0;i<splitN;i++) totHistogram[i]=0.0;

      boost::mpi::all_reduce(*mpiWorld, histogram, splitN, totHistogram, std::plus<real>());
      boost::mpi::all_reduce(*mpiWorld, num_part, total_num_part, std::plus<int>());

      // normalizing
      int nconfigs = 1; //config - 1
      real rho = (real)total_num_part / splitN;
      
      for(int i=0; i < splitN; i++){
        totHistogram[i] /= rho;
      }

      python::list pyli;
      for(int i=0; i < splitN; i++){
        pyli.append( totHistogram[i] );
      }
      // totHistogram.clear();
      // histogram.clear();
      delete[] totHistogram;
      totHistogram = 0;
      delete[] histogram;
      histogram = 0;

      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real XDensity::compute() const {
      return -1.0;
    }
    
    using namespace boost::python;

    void XDensity::registerPython() {
      using namespace espressopp::python;
      class_<XDensity, bases< Observable > >
        ("analysis_XDensity", init< shared_ptr< System > >())
        .def("compute", &XDensity::computeArray)
      ;
    }
  }
}
