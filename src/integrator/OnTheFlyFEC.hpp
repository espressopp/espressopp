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
#ifndef _INTEGRATOR_ONTHEFLYFEC_HPP
#define _INTEGRATOR_ONTHEFLYFEC_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Real3D.hpp"
#include "SystemAccess.hpp"
#include "interaction/Interpolation.hpp"
#include <map>

#include "python.hpp"
#include "types.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"


namespace espressopp {
  namespace integrator {

    class OnTheFlyFEC : public Extension {

      public:
        OnTheFlyFEC(shared_ptr<System> system);
        virtual ~OnTheFlyFEC();

        void setBins(int bins);
        int getBins();
        
        void setGap(int gap);
        int getGap();
        
        void setSteps(int steps);
        int getSteps();
        
        python::list writeFEC();
        
        void makeArrays();
        
        void resetCounter();

        void setCenter(real x, real y, real z);
        
        static void registerPython();

      private:

        boost::signals2::connection _gatherStats;

        int bins;
        int gap;
        int steps;
   
        int counter;
        int gapcounter;
        
        //std::vector<int> NumbersAtoms;
        //std::vector<real> EnergyDiff;
        real * EnergyDiff;
        int * NumbersAtoms;
        
        void gatherStats();
        
        void connect();
        void disconnect();

        Real3D center;
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif