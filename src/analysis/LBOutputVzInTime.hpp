/*
 Copyright (C) 2012-2016
     Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011
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

#ifndef _LBOUTPUT_VZINTIME_HPP
#define _LBOUTPUT_VZINTIME_HPP

#include "LBOutput.hpp"

namespace espressopp {
   namespace analysis {
      class LBOutputVzInTime : public LBOutput {
      public:
         LBOutputVzInTime(shared_ptr<System> _system,
                          shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann);

         void writeOutput();

         static void registerPython();
      };
   }
}

#endif
