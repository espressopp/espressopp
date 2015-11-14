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
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espressopp {

  namespace integrator {

      class Adress : public Extension {

      public:
        shared_ptr<VerletListAdress> verletList;
        shared_ptr<FixedTupleListAdress> fixedtupleList;
        bool KTI;
        
        real dhy;
        real pidhy2;
        real dex;
        real dex2;
        real dexdhy;
        real dexdhy2;
        
        Adress(shared_ptr<System> _system, shared_ptr<VerletListAdress> _verletList, shared_ptr<FixedTupleListAdress> _fixedtupleList, bool _KTI = false);

        ~Adress();
        
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        boost::signals2::connection _SetPosVel, _initForces, _integrate1, _inIntP, _integrate2, _recalc2, _befIntV;  //_aftCalcF;
        
        void integrate1(real&);
        void initForces();
        void SetPosVel();
        void integrate2();
        void aftCalcF();
        void communicateAdrPositions();

        void connect();
        void disconnect();
        
        real weight(real);
        real weightderivative(real);
      };

  }

}

#endif
