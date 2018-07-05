/*
  Copyright (C) 2012-2018
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

// ESPP_CLASS
#ifndef _INTEGRATOR_FREEENERGYCOMPENSATION_HPP
#define _INTEGRATOR_FREEENERGYCOMPENSATION_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Real3D.hpp"
#include "SystemAccess.hpp"
#include "interaction/Interpolation.hpp"
#include <unordered_map>


#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"


namespace espressopp {
  namespace integrator {

    class FreeEnergyCompensation : public Extension {

      public:
        bool sphereAdr;
        int ntrotter;
        bool slow;
        FreeEnergyCompensation(shared_ptr<System> system, bool _sphereAdr = false, int _ntrotter = 1, bool _slow = false);
        virtual ~FreeEnergyCompensation();

        /** Setter for the filename, will read in the table. */
        void addForce(int itype, const char* _filename, int type);
        const char* getFilename() const { return filename.c_str(); }

        void applyForce();

        real computeCompEnergy();

        void setCenter(real x, real y, real z);

        static void registerPython();

      private:

        boost::signals2::connection _applyForce;

        void connect();
        void disconnect();

        Real3D center;
        std::string filename;
        typedef shared_ptr <interaction::Interpolation> Table;
        std::unordered_map<int, Table> forces; // map type to force

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
