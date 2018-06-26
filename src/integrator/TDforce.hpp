/*
  Copyright (C) 2017,2018
      Jakub Krajniak (jkrajniak at gmail.com), Max Planck Institute for Polymer Research
  Copyright (C) 2012,2013,2014,2015,2016
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
#ifndef _INTEGRATOR_TDFORCE_HPP
#define _INTEGRATOR_TDFORCE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Real3D.hpp"
#include "SystemAccess.hpp"
#include "VerletListAdress.hpp"
#include "interaction/Interpolation.hpp"
#include <unordered_map>


#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"


namespace espressopp {
  namespace integrator {

    class TDforce : public Extension {

      public:
        shared_ptr<VerletListAdress> verletList;
        real startdist;
        real enddist;
        int edgeweightmultiplier;
        bool slow;
        TDforce(shared_ptr<System> system, shared_ptr<VerletListAdress> _verletList, real _startdist = 0.0, real _enddist = 0.0, int _edgeweightmultiplier = 1, bool _slow = false);

        ~TDforce();

        /** Setter for the filename, will read in the table. */
        void addForce(int itype, const char* _filename, int type);
        const char* getFilename() const { return filename.c_str(); }

        void applyForce();
        real computeTDEnergy();

        static void registerPython();

      private:

        boost::signals2::connection _applyForce;

        void connect();
        void disconnect();
        real getForce(longint type_id, real dist);

        Real3D center; // center of adress zone, from verletlistadress (assumes only one point as center)
        bool sphereAdr; // true: adress region is spherical centered on point x,y,z or particle pid; false: adress region is slab centered on point x or particle pid, from verletlistadres
        std::string filename;
        typedef shared_ptr <interaction::Interpolation> Table;
        std::unordered_map<int, Table> forces; // map type to force

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
