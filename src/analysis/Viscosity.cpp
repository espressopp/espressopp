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

#include "bc/BC.hpp"
#include "Viscosity.hpp"
#include "esutil/Error.hpp"
#include "PressureTensor.hpp"

using namespace std;

namespace espressopp {
  namespace analysis {

    void Viscosity::gather() {
      Tensor prt = PressureTensor(getSystem()).computeRaw();
      Autocorrelation::gather( Real3D(prt[3], prt[4], prt[5]) );
    }
    
    python::list Viscosity::compute(real t0, real dt, real T) {
      python::list auto_pxy_pxy_py = Autocorrelation::compute();
      size_t M = getListSize();
      
      cout << "Size of array: " << M << endl;
      
      real* auto_pxy_pxy = new real[M];
      
      System& system = getSystemRef();
      int n_nodes = system.comm->size();
      int this_node = system.comm->rank();
      
      Real3D Li = system.bc->getBoxL();
      // factor: V/T * integral <Pxy * Pxy>
      real V_T = Li[0] * Li[1] * Li[2] / T;
      
      if(this_node==0){
        for (int i = 0; i < len(auto_pxy_pxy_py); ++i){
          auto_pxy_pxy_py[i] = V_T * boost::python::extract<double>(auto_pxy_pxy_py[i]);
        }
      }
      
      boost::mpi::broadcast(*system.comm, auto_pxy_pxy, (int)M, 0);
    
      // load distribution
      // it is 1+2+3+...+M
      real local_num = ( M * (M + 1) / 2.0) / (real)(n_nodes) + 1.0;
      real num_m [n_nodes];
      for(int i=0; i<n_nodes; i++){
        num_m[i] = 0.0;
        real max_num = (i+1) * local_num;
        int m_max = int( ( sqrt(1.0+8.0*max_num) - 1.0 )/2. );
        int lastNum = 0;
        if(i!=0)
          lastNum = num_m[i-1];
                  
        if(m_max-lastNum==0)
          num_m[i] = num_m[i-1] + 1;
        else
          num_m[i] = m_max;

        if(num_m[i] >= M)
          num_m[i] = M - 1;
      }
            
      int min_m = 1;
      if(this_node!=0)
        min_m = num_m[this_node-1];

      int max_m = num_m[this_node];

      python::list integr;

      if(this_node==0)
        cout<<"integration..." <<endl;

      real SUM = (auto_pxy_pxy[0]+auto_pxy_pxy[min_m]) * 0.5;
      for (int j=1; j<min_m; j++)
        SUM += auto_pxy_pxy[j];

      integr.append( python::make_tuple(t0+min_m*dt, SUM*dt) );
      for(int i=min_m+1; i<max_m; i++){
        SUM += (auto_pxy_pxy[i-1]+auto_pxy_pxy[i]) * 0.5;
        integr.append( python::make_tuple(t0+i*dt, SUM*dt) );
      }

      delete [] auto_pxy_pxy;
      auto_pxy_pxy = NULL;
      return integr;
    }
    
    // Python wrapping
    void Viscosity::registerPython() {
      using namespace espressopp::python;

      class_<Viscosity>( "analysis_Viscosity", init< shared_ptr< System > >() )
      .def("gather", &Viscosity::gather)
      .def("compute", &Viscosity::compute)
      ;
    }
  }
}
