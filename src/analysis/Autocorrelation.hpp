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
#ifndef _ANALYSIS_AUTOCORRELATION_HPP
#define _ANALYSIS_AUTOCORRELATION_HPP

#include "python.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

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

    class Autocorrelation : public SystemAccess {

    public:
      // Constructor, allow for unlimited snapshots.
      Autocorrelation(shared_ptr<System> system): SystemAccess (system){
      }
      ~Autocorrelation(){
        valueList.clear();
      }

      // get number of available snapshots. Returns the size of ValueList
      unsigned int getListSize() const;

      // Take a snapshot (save the current value)
      void gather(Real3D);
      
      // Get a configuration from ConfigurationList
      Real3D getValue(unsigned int position) const;

      // it returns all values
      vector<Real3D> all() const;

      // it erases all the configurations from ConfigurationList
      void clear(){
        valueList.clear();
      }

      python::list compute();

      static void registerPython();
    
    private:

      void pushValue(Real3D);
 
      // the list of snapshots
      vector<Real3D> valueList;
    };
  }
}

#endif
