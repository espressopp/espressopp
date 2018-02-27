/*
  Copyright (C) 2012,2013,2017,2018
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

#include "python.hpp"
#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "mpi.hpp"
#include "Particle.hpp"
#include "../../Plumed.h"

namespace espressopp {
  namespace integrator {

    /** ExtPlumed */
    class ExtPlumed : public Extension {

    public:
      ExtPlumed(shared_ptr < System >, python::object, std::string, std::string, real);
      void applyForceToAll();
      void updatePlumed();
      real getBias();
      void setUnitStyle(std::string);
      void setTimeUnit(real);
      void setEnergyUnit(real);
      void setLengthUnit(real);
      void Init();
      bool getChargeState();
      void setChargeState(bool);
      virtual ~ExtPlumed();
      /** Register this class so it can be used from Python. */
      static void registerPython();

    private:
      // pointer to plumed object:
      PLMD::Plumed * p;
      std::string plumedfile;
      std::string units;
      std::string plumedlog;
      real dt;

      longint nreal; // total number of atoms (real & ghost) on the processor
      longint natoms; // total number of atoms
      int *  gatindex;
      real * masses;
      real * charges;
      real * pos, *f;

      boost::signals2::connection _aftCalcF;
      boost::signals2::connection _aftIntV;

      void connect();
      void disconnect();
      // void getTimeStep();
      /** Logger */
      static LOG4ESPP_DECL_LOGGER(theLogger);
      bool chargeState;
      real bias;
    };
  }
}

#endif
