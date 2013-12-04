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
#ifndef _ADRESS_HPP
#define _ADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espresso {

  namespace integrator {

      class Adress : public Extension {

      public:

        Adress(shared_ptr<System> system);

        ~Adress();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        boost::signals2::connection _initForces, _integrate1, _integrate2;

        void integrate1(real&);
        void initForces();
        void integrate2();

        void connect();
        void disconnect();
      };

  }

}

#endif
