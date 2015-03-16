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
#ifndef _INTEGRATOR_EXTANALYZE_HPP
#define _INTEGRATOR_EXTANALYZE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
//#include "analysis/AnalysisBase.hpp"
#include "ParticleAccess.hpp"

namespace espressopp {
  //using namespace analysis;
  namespace integrator {

    /** ExtAnalyze */
    class ExtAnalyze : public Extension {
      public:
        //ExtAnalyze(shared_ptr< AnalysisBase > _analysis, int _interval);
        ExtAnalyze(shared_ptr< ParticleAccess > _particle_access, int _interval);
        virtual ~ExtAnalyze() {};
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        boost::signals2::connection _aftIntV;
        void connect();
        void disconnect();
        void perform_action();
        //void performMeasurement();

        shared_ptr< ParticleAccess > particle_access;
        int interval;
        int counter;

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
