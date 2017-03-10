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

#include "python.hpp"
#include "CMVelocity.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "mpi.h"

using namespace espressopp;

namespace espressopp {
  namespace analysis {

    using namespace iterator;

    LOG4ESPP_LOGGER(CMVelocity::logger, "CMVelocity");

    Real3D CMVelocity::computeRaw() {

      System& system = getSystemRef();
  
      real myMass;
      Real3D myV;

      CellList realCells = system.storage->getRealCells();

      // Compute the local velocity
      myMass=0.;
      myV=0.;
      
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D& velocity = cit->velocity();
        real mass = cit->getMass();
        myMass += mass;
        myV += mass*velocity;
      }

      real v_data[4], total_v[4];
      v_data[0] = myMass;
      v_data[1] = myV[0];
      v_data[2] = myV[1];
      v_data[3] = myV[2];

      // Compute the total velocity of the system
      boost::mpi::all_reduce(*(system.comm), v_data, 4, total_v, std::plus<real>());

      v[0] = total_v[1]/total_v[0];
      v[1] = total_v[2]/total_v[0];
      v[2] = total_v[3]/total_v[0];

      return v;
    }

    void CMVelocity::reset() {
      System& system = getSystemRef();

      computeRaw();
      CellList realCells = system.storage->getRealCells();
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        cit->velocity() -= v;
      }
      
      // set private variable to zero (new CM-velocity).
      // one could use computeRaw() again, but it is slower
      v = Real3D(0.);

    }

    Real3D CMVelocity::getV() const { return v; }
    
    // Python wrapping

    void CMVelocity::registerPython() {

      using namespace espressopp::python;

      class_<CMVelocity, bases<AnalysisBase> >
        ("analysis_CMVelocity", init<shared_ptr<System> >())
      .add_property("v", &CMVelocity::getV)
      ;
    }
  }
}
