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
#ifndef _ANALYSIS_INTRACHAINDISTSQ_HPP
#define _ANALYSIS_INTRACHAINDISTSQ_HPP

#include "python.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"
#include "FixedPairList.hpp"
#include "AllParticlePos.hpp"

namespace espressopp {
  namespace analysis {

    /** calculate mean squared intra-chain-distances
    */

    class IntraChainDistSq : public AllParticlePos {
    public:
      shared_ptr<FixedPairList> fpl;
      IntraChainDistSq(shared_ptr<System> system, shared_ptr<FixedPairList> _fpl) : AllParticlePos(system), fpl(_fpl) {};
      ~IntraChainDistSq() {};
      python::list compute();
      static void registerPython();

    protected:
      static LOG4ESPP_DECL_LOGGER(logger);
    };
  }
}

#endif
