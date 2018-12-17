/*
  Copyright (C) 2012,2013,2018
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
#include "AngularPotential.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(AngularPotential::theLogger, "AngularPotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void AngularPotential::registerPython() {
        using namespace espressopp::python;

        real (AngularPotential::*computeEnergy1)
            (const Real3D& dist12, const Real3D& dist32) const =
                &AngularPotential::computeEnergy;

        real (AngularPotential::*computeEnergy2)
            (real theta) const =
                &AngularPotential::computeEnergy;

        void (AngularPotential::*computeForce1)
            (Real3D& force12, Real3D& force32,
            const Real3D& dist12, const Real3D& dist32) const =
                &AngularPotential::computeForce;

        real (AngularPotential::*computeForce2)(real theta) const =
                &AngularPotential::computeForce;

        class_< AngularPotential, boost::noncopyable >
            ("interaction_AngularPotential", no_init)
            .add_property("cutoff",
                &AngularPotential::getCutoff,
                &AngularPotential::setCutoff)
            .add_property("colVarBondList",
                &AngularPotential::getColVarBondList,
                &AngularPotential::setColVarBondList)
            .add_property("colVarAngleList",
                          &AngularPotential::getColVarAngleList,
                          &AngularPotential::setColVarAngleList)
            .add_property("colVarDihedList",
                &AngularPotential::getColVarDihedList,
                &AngularPotential::setColVarDihedList)
            .add_property("colVar",
                &AngularPotential::getColVar)
            .add_property("addCVBond",
                &AngularPotential::getColVarBondList,
                &AngularPotential::setColVarBondList)
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce1))
            .def("computeForce", pure_virtual(computeForce2))
        ;
    }
  }
}
