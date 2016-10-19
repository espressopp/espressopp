/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research
  
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
#ifndef _RATTLE_HPP
#define _RATTLE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include "boost/signals2.hpp"

namespace espressopp {
  namespace integrator {
    class Rattle : public Extension {

      public:
        Rattle(shared_ptr<System> _system, real _maxit, real _tol, real _rptol);
        ~Rattle();        

        void addBond(int pid1, int pid2, real constraintDist, real mass1, real mass2);

        void saveOldPos();
        void applyPositionConstraints();
        void applyVelocityConstraints();

        static void registerPython();

      private:
        boost::signals2::connection _befIntP, _aftIntP, _aftIntV;
        void connect();
        void disconnect();


        // positions in previous timestep
        // key is pid of light atom in bond, value is heavy atom and light atom coordinates
        typedef boost::unordered_map<longint, std::pair<Real3D, Real3D> > OldPos;
        OldPos oldPos;

        std::vector<longint> lightPart; //list of light particles in constrained bonds on the CPU, for use in applyVelocityConstraints

        struct ConstrainedBond {
          longint pidHeavy;
          longint pidHyd;
          real constraintDist2; //squared distance
          real invmassHeavy; 
          real invmassHyd;
        };
        boost::unordered_map<longint, ConstrainedBond> constrainedBonds; //key: light atom pid
        std::vector<longint> constrainedBondsKeys; //list of keys (light atom pids) in constrainedBonds map

        real maxit; //maximum number of iterations
        real tol; //tolerance for deciding if constraint distance and current distance are similar enough
        real rptol; //tolerance for deciding if the angle between the bond vector at end of previous timestep and current vector has become too large
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
