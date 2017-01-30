/*
 Copyright (C) 2012-2016 Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

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
#ifndef _ANALYSIS_LBOUTPUT_HPP
#define _ANALYSIS_LBOUTPUT_HPP

#include "integrator/LatticeBoltzmann.hpp"
#include "AnalysisBase.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
   namespace analysis {
      using namespace iterator;
      /** Abstract base class for arbitrary output from LB simulations. */
      class LBOutput : public AnalysisBaseTemplate< real > {
      public:
         /* Constructor for the class */
         LBOutput(shared_ptr< System > _system,
                  shared_ptr< integrator::LatticeBoltzmann > _latticeboltzmann)
         : AnalysisBaseTemplate< real >(_system) {
            latticeboltzmann = _latticeboltzmann;
         }
         /* Destructor for the class */
         virtual ~LBOutput () {}

         real computeRaw() {
            writeOutput();
            real _dummy = 0.;
            return _dummy;
         }

         python::list compute() {
            python::list _dummy;
            real res = computeRaw();
            _dummy.append(0.);
            return _dummy;
         }

         python::list getAverageValue() {
            python::list _dummy;
            _dummy.append(0.);
            return _dummy;
         }

         void resetAverage() {}
         void updateAverage(real res) {return;}

         /* writing of a profile into the output */
         virtual void writeOutput () = 0;

         static void registerPython();
      protected:
         shared_ptr<integrator::LatticeBoltzmann> latticeboltzmann;
      };
   }
}

#endif
