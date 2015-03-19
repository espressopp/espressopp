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
#include "RadialDistrF.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"

#include <boost/serialization/map.hpp>

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {
    
    // TODO currently works correctly for Lx = Ly = Lz
    // rdfN is a level of discretisation of rdf (how many elements it contains)
    python::list RadialDistrF::computeArray(int rdfN) const {

      System& system = getSystemRef();
      esutil::Error err(system.comm);
      Real3D Li = system.bc->getBoxL();
      Real3D Li_half = Li / 2.;
      
      int nprocs = system.comm->size();
      int myrank = system.comm->rank();
      
      real histogram [rdfN];
      for(int i=0;i<rdfN;i++) histogram[i]=0;
          
      real dr = Li_half[1] / (real)rdfN; // If you work with nonuniform Lx, Ly, Lz, you
                                         // should use for Li_half[XXX] the shortest side length   
      
      int num_part = 0;
      ConfigurationPtr config = make_shared<Configuration> ();
      for (int rank_i=0; rank_i<nprocs; rank_i++) {
        map< size_t, Real3D > conf;
        if (rank_i == myrank) {
          CellList realCells = system.storage->getRealCells();
          for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            int id = cit->id();
            conf[id] = cit->position();
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
      }
      // now all CPUs have all particle coords and num_part is the total number of particles
      
      // use all cpus
      // TODO it could be a problem if   n_nodes > num_part
      int numi = num_part / nprocs + 1;
      int mini = myrank * numi;
      int maxi = mini + numi;
      
      if(maxi>num_part) maxi = num_part;
      
      int perc=0;
      real denom = 100.0 / (real)(maxi-mini);
      for(int i = mini; i<maxi; i++){
        Real3D coordP1 = config->getCoordinates(i);
        for(int j = i+1; j<num_part; j++){
          Real3D coordP2 = config->getCoordinates(j);
          Real3D distVector = coordP1 - coordP2;
          
          // minimize the distance in simulation box
          for(int ii=0; ii<3; ii++){
            if( distVector[ii] < -Li_half[ii] ) distVector[ii] += Li[ii];
            if( distVector[ii] >  Li_half[ii] ) distVector[ii] -= Li[ii];
          }
          
          int bin = (int)( distVector.abs() / dr );
          if( bin < rdfN){
            histogram[bin] += 1.0;
          }
          
        }
        /*
         * additional calculations slow down routine but from the other hand
         * it helps to monitor progress
         */
        if(print_progress && system.comm->rank()==0){
          perc = (int)((i-mini)*denom);
          if(perc%5==0){
            cout<<"calculation progress (radial distr. func.): "<< perc << " %\r"<<flush;
          }
        }
        
      }
      if(system.comm->rank()==0)
        cout<<"calculation progress (radial distr. func.): 100 %" <<endl;

      real totHistogram[rdfN];
      boost::mpi::all_reduce(*system.comm, histogram, rdfN, totHistogram, plus<real>());
      
      // normalizing
      int nconfigs = 1; //config - 1
      //avg_part = part_total / float(nconfigs)
      real rho = (real)num_part / (Li[0]*Li[1]*Li[2]);
      real factor = 2.0 * M_PIl * dr * rho * (real)nconfigs * (real)num_part;
      
      for(int i=0; i < rdfN; i++){
        real radius = (i + 0.5) * dr;
        totHistogram[i] /= factor * (radius*radius + dr*dr / 12.0);
      }
      
      python::list pyli;
      for(int i=0; i < rdfN; i++){
        pyli.append( totHistogram[i] );
      }
      return pyli;
    }

    // TODO: this dummy routine is still needed as we have not yet ObservableVector
    real RadialDistrF::compute() const {
      return -1.0;
    }

    void RadialDistrF::registerPython() {
      using namespace espressopp::python;
      class_<RadialDistrF, bases< Observable > >
        ("analysis_RadialDistrF", init< shared_ptr< System > >())
        .add_property("print_progress", &RadialDistrF::getPrint_progress, &RadialDistrF::setPrint_progress)
        .def("compute", &RadialDistrF::computeArray)
      ;
    }
  }
}
