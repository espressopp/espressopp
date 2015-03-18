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
#ifndef _ANALYSIS_NEIGHBORFLUCTUATION_HPP
#define _ANALYSIS_NEIGHBORFLUCTUATION_HPP

#include "python.hpp"
#include "types.hpp"
#include "Observable.hpp"
#include "esutil/Error.hpp"
#include "bc/BC.hpp"
#include "Real3D.hpp"

#include "storage/DomainDecomposition.hpp"
#include "storage/CellGrid.hpp"

namespace espressopp {
  namespace analysis {
    /** Class to get the number of particles in the system. */
    class NeighborFluctuation : public Observable {
    public:
      NeighborFluctuation(shared_ptr< System > system, real _radius) : Observable(system), radius(_radius){
        esutil::Error err(system->comm);
        
        Real3D Li = system->bc->getBoxL();
        
        Int3D cellG = system->storage->getInt3DCellGrid();
        
        real minL = std::min(Li[0]/(real)cellG[0], std::min(Li[1]/(real)cellG[1],Li[2]/(real)cellG[2]));
        if(radius > minL){
          std::stringstream msg;
          msg<<"Error. Radius for checking near neighbors should be smaller then "
                  "minimal cell size.\n";
          msg<<"Otherwise it might count the same particle twice or miss some particles. "
                  "radius="<<radius<<" minCellSize="<<minL;
          err.setException( msg.str() );
        }
      }
      virtual ~NeighborFluctuation() {}
      virtual real compute() const;
      
      // it returns two values <n^2>-<n>^2 and <n>
      virtual python::list computeValues() const;

      static void registerPython();

    protected:
      real radius;
    };
  }
}

#endif
