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
#ifndef _ANALYSIS_ALLPARTICLEPOS_HPP
#define _ANALYSIS_ALLPARTICLEPOS_HPP

#include "types.hpp"
#include "SystemAccess.hpp"
#include <map>

namespace espressopp {
  namespace analysis {

    /** Class provides pid and position information of all particles on all cpus
        (used e.g. for atom decomposition parallel analysis)
    */

    struct sBuf {
      real r[3];
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
   	    for (int i = 0; i < 3; ++i) ar & r[i];
      }
    };

    typedef std::map< size_t, sBuf > ConfMap;

    class AllParticlePos : public SystemAccess {
    public:
      AllParticlePos(shared_ptr<System> system) : SystemAccess (system) {};
      ~AllParticlePos() {};
      /** gather and broadcast all particle positions to all cpus */
      void gatherAllPositions();
      
      ConfMap getAllPositions();
      
      static void registerPython();

      ConfMap AllPositions;
      int numParticles;

    protected:
      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}

#endif
