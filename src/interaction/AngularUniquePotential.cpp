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

#include "python.hpp"
#include "AngularUniquePotential.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(AngularUniquePotential::theLogger, "AngularUniquePotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void AngularUniquePotential::registerPython() {
      using namespace espressopp::python;

      real (AngularUniquePotential::*computeEnergy1)
          (const Real3D& r12, const Real3D& r32, real theta0) const =
              &AngularUniquePotential::computeEnergy;

      real (AngularUniquePotential::*computeEnergy2)
          (real theta, real theta0) const =
              &AngularUniquePotential::computeEnergy;

      void (AngularUniquePotential::*computeForce1)
          (Real3D& force12, Real3D& force32,
          const Real3D& r12, const Real3D& r32, real theta0) const =
              &AngularUniquePotential::computeForce;

      real (AngularUniquePotential::*computeForce2)(real theta, real theta0) const =
              &AngularUniquePotential::computeForce;

      class_< AngularUniquePotential, boost::noncopyable >
          ("interaction_AngularUniquePotential", no_init)
          .add_property("cutoff",
              &AngularUniquePotential::getCutoff,
              &AngularUniquePotential::setCutoff)
          .def("computeEnergy", pure_virtual(computeEnergy1))
          .def("computeEnergy", pure_virtual(computeEnergy2))
          .def("computeForce", pure_virtual(computeForce1))
          .def("computeForce", pure_virtual(computeForce2))
      ;
    }
  }
}
