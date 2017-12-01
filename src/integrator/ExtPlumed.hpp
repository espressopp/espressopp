/*
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
#ifndef _INTEGRATOR_ExtPlumed_HPP
#define _INTEGRATOR_ExtPlumed_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"

#include "../../Plumed.h"

namespace espressopp {
  namespace integrator {

    /** ExtPlumed */

    class ExtPlumed : public Extension {

      public:
      ExtPlumed(shared_ptr < System >, std::string, std::string, std::string);
      void applyForceToAll();
        real getBias();
        virtual ~ExtPlumed() {};
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        // pointer to plumed object:
        PLMD::Plumed*p;
        std::string plumedfile;
        std::string units;
      std::string plumedlog;
      // real peTotal;
        real dt;


        longint nlocal; // total number of atoms (real & ghost) on the processor
        longint natoms; // total number of atoms
        size_t *gatindex;
        real *masses;
        real *forces;
        real *pos;
        real *charges;

        boost::signals2::connection _aftCalcF;
        void connect();
        void disconnect();
        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
        const bool charged;
        real bias;
    };
  }
}

#endif
