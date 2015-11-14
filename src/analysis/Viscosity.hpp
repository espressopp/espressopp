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

// ESPP_CLASS
#ifndef _ANALYSIS_VISCOSITY_HPP
#define _ANALYSIS_VISCOSITY_HPP

#include "types.hpp"
#include "Autocorrelation.hpp"
#include "python.hpp"

using namespace std;

namespace espressopp {
  namespace analysis {

    /*
     * Class stores some single value in time for later calculation (using parallel
     * computation) of autocorrelation function. It is useful for example for viscosity
     * calculations.
     * 
     * !Important! It should be the same time period between snapshots.
    */
    
    // now the single value is Real3D
    // TODO probably template realization

    class Viscosity : public Autocorrelation {

    public:
      // Constructor, allow for unlimited snapshots.
      Viscosity(shared_ptr<System> system): Autocorrelation (system){
      }
      ~Viscosity(){
      }

      // Take a snapshot (save the current value of nonlinar component of pressure tensor)
      void gather();
      
      python::list compute(real t0, real dt, real T);

      static void registerPython();
      
    };
  }
}

#endif
