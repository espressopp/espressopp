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

#include "Autocorrelation.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/Error.hpp"
#include "mpi.h"

using namespace std;

namespace espressopp {
  namespace analysis {

    using namespace iterator;

    unsigned int Autocorrelation::getListSize() const{
      return valueList.size();
    }
   
    vector<Real3D> Autocorrelation::all() const{
      return valueList;
    }

    Real3D Autocorrelation::getValue(unsigned int position) const{
      unsigned int nconfigs = getListSize();
      if (0 <= position and position < nconfigs) {
        return valueList[position];
      }
      else{
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        stringstream msg;
        msg << "Error. Velocities::get <out-of-range>";
        err.setException( msg.str() );
        return Real3D(0,0,0);
      }
    }

    void Autocorrelation::pushValue(Real3D value){
      valueList.push_back(value);
    }

    void Autocorrelation::gather(Real3D value) {
      pushValue(value);
    }
    
    python::list Autocorrelation::compute() {
      size_t M = getListSize();
      
      System& system = getSystemRef();
      
      int n_nodes = system.comm -> size();
      int this_node = system.comm -> rank();

      // TODO it could be a problem if   n_nodes > total_num !!!
      unsigned int num_m[n_nodes];
      unsigned int num_mH[n_nodes];
      
      if(this_node == 0){
        
        // it is 1+2+3+...+M
        double local_num = ( (double)M * (double)(M + 1) / 2.0) / (double)n_nodes + 1.0;
        
        for(int i=0; i<n_nodes; i++){
          double max_num = (i+1) * local_num;
          
          unsigned int m_max = (unsigned int) ( ( sqrt(1.0+8.0*max_num) - 1.0 )/2. );
          
          unsigned int lastNum = (i==0) ? 0 : num_mH[i-1];
          
          if(m_max-lastNum==0)
            num_mH[i] = num_mH[i-1] + 1;
          else
            num_mH[i] = m_max;
          
          if(num_mH[i] > M) num_mH[i] = M;
        }
        
        for(int i=0; i<n_nodes-1; i++)
          num_m[i] = M - num_mH[n_nodes-2-i];
        num_m[n_nodes-1] = M;
      }

      boost::mpi::broadcast(*system.comm, num_m, n_nodes, 0);
      
      // now num_m[i], i - cpu number, is a number of series in "for" statement for each cpu
      
      unsigned int min_m = (this_node==0) ? 0 : num_m[this_node-1];
      unsigned int max_m = num_m[this_node];

      real* Z;
      Z = new real[M];
      
      cout<< "calculating autocorrelation.." << endl;
      int perc=0;
      real denom = 100.0 / (real)(max_m-min_m);
      for(unsigned int m=min_m; m<max_m; m++){
        Z[m] = 0.0;
        for(unsigned int n=0; n<M-m; n++){
          Z[m] += getValue(n + m) *  getValue(n);
        }
        
        /*
         * additional calculations slow down routine but from the other hand
         * it helps to monitor progress
         */
        if(system.comm->rank()==0){
          perc = (int)(m*denom);
          if(perc%5==0){
            cout<<"calculation progress (autocorrelation): "<< perc << " %\r" <<flush;
          }
        }
      }
      
      if(system.comm->rank()==0)
        cout<<"calculation progress (autocorrelation): 100 %" <<endl;
      
      real coef = 3.0; // only if value is Real3D
      
      for(unsigned int m=min_m; m<max_m; m++){
        Z[m] /= ( (real)(M-m)*coef );
      }
      
      // TODO probably could be done nicer. gather doesn't work with different length of array
      unsigned long int MM = M * n_nodes;
      real* totZ = new real[MM];
      boost::mpi::gather(*system.comm, Z, M, totZ, 0);
      
      python::list pyli;
      
      if(this_node == 0){
        int count = 0;
        for(int m=0; m<M; m++){
          if( m >= num_m[count] ) count ++;
          pyli.append( totZ[M*count + m] );
        }
      }
      
      delete [] Z;
      Z = NULL;
      delete [] totZ;
      totZ = NULL;
      
      return pyli;
    }
    
    // Python wrapping
    void Autocorrelation::registerPython() {
      using namespace espressopp::python;

      class_<Autocorrelation>(
        "analysis_Autocorrelation",
        init< shared_ptr< System > >()
      )
      .def_readonly("size", &Autocorrelation::getListSize)
      
      .def("gather", &Autocorrelation::gather)
      .def("__getitem__", &Autocorrelation::getValue)
      .def("all", &Autocorrelation::all)
      .def("clear", &Autocorrelation::clear)
      .def("compute", &Autocorrelation::compute)
      
      ;
    }
  }
}
