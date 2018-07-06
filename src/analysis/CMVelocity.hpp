/*
  Copyright (C) 2014 Pierre de Buyl
  Copyright (C) 2012,2013,2017
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
#ifndef _ANALYSIS_TOTAL_VELOCITY_HPP
#define _ANALYSIS_TOTAL_VELOCITY_HPP

#include "types.hpp"
#include "Real3D.hpp"
#include "AnalysisBase.hpp"

namespace espressopp {
  namespace analysis {

    /** Class to compute the total velocity of a system. CMVelocity provides
	a facility to reset the total velocity of the system.
    */

    class CMVelocity : public AnalysisBaseTemplate <Real3D> {

    public:

      CMVelocity(shared_ptr<System> system) : AnalysisBaseTemplate <Real3D> (system) {}

      ~CMVelocity() {}

      /** Reset the total velocity of the system*/
      void reset();

      void perform_action() { reset(); }

      /** get current velocity */
      Real3D getV() const;

      static void registerPython();

      /** virtual functions from AnalysisBase class */
      Real3D computeRaw(); // compute the total velocity of the system

      python::list compute() {
        python::list ret;
        Real3D res = computeRaw();
        ret.append(res);
        return ret;
      }
      
      python::list getAverageValue() {
        python::list ret;
        return ret;
      }
      
      void resetAverage() {return;}

      void updateAverage(Real3D res) {return;}
      
    protected:

      static LOG4ESPP_DECL_LOGGER(logger);

    private:

      Real3D v;
    };
  }
}

#endif
