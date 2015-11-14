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
#ifndef _ANALYSIS_TEST_HPP
#define _ANALYSIS_TEST_HPP

#include "types.hpp"
#include "AnalysisBase.hpp"

namespace espressopp {
  namespace analysis {
    /** Class to test AnalysisBase */
    class Test : public AnalysisBaseTemplate< int > {
    public:
      static void registerPython();

      Test(shared_ptr< System > system) : AnalysisBaseTemplate< int >(system) {};
      virtual ~Test() {};

      int computeRaw() {
        return 99;
      }

      python::list compute() {
        python::list ret;
        int res = computeRaw();
        ret.append(res);
        return ret;
      }

      python::list getAverageValue() {
        python::list ret;
        int res;
        res = nMeasurements>0 ? newAverage : 0;
        ret.append(res);
        res = nMeasurements>0 ? newVariance : 0;
        ret.append(res);
        return ret;
      }

      void resetAverage() {
        newAverage = 0.0;
        newVariance = 0.0;
      }

      void updateAverage(int res) {
        // compare e. g. Knuth TAOCP vol 2, 3rd edition, page 232
        if (nMeasurements == 1) {
            newAverage     = res;
            lastAverage    = newAverage;
        } else {
            newAverage  = lastAverage  + (res - lastAverage) / nMeasurements;
            newVariance = lastVariance + (res - lastAverage) * (res - newAverage);
            lastAverage = newAverage;
            lastVariance = newVariance;
        }
        return;
      }
    };
  }
}

#endif
