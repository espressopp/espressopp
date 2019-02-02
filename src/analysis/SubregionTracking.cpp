/*
  Copyright (C) 2017,2018
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

#include "python.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "SubregionTracking.hpp"
#include "bc/BC.hpp"

#include <boost/serialization/map.hpp>

using namespace espressopp;
using namespace iterator;
using namespace std;

namespace espressopp {
  namespace analysis {
    int SubregionTracking::compute_int() const {

      int myN, systemN;
      myN = 0;
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        if(particlelist.count(cit->id())) {
          Particle &vp = *cit;
          Real3D dist(0.0, 0.0, 0.0);
          system.bc->getMinimumImageVector(dist, vp.position(), center);

          if(geometry == spherical) {
            if (dist.abs() <= span) myN += 1;
          }
          else if(geometry == x_bounded) {
            if (fabs(dist[0]) <= span) myN += 1;
          }
          else if(geometry == y_bounded) {
            if (fabs(dist[1]) <= span) myN += 1;
          }
          else if(geometry == z_bounded) {
            if (fabs(dist[2]) <= span) myN += 1;
          }
          else {
            throw std::invalid_argument("Invalid geometry parameter.");
          }
        }
      }

      boost::mpi::reduce(*getSystem()->comm, myN, systemN, std::plus<int>(), 0);
      return systemN;
    }

    void SubregionTracking::setCenter(real x, real y, real z) {
      center = Real3D(x, y, z);
    }

    void SubregionTracking::registerPython() {
      using namespace espressopp::python;
      void (SubregionTracking::*pySetCenter)(real x, real y, real z) = &SubregionTracking::setCenter;
      class_<SubregionTracking, bases< Observable > >
      ("analysis_SubregionTracking", init< shared_ptr< System >, real, int >())
      .def("setCenter", pySetCenter)
      .def("addPID", &SubregionTracking::addPID)
      ;
    }
  }
}
