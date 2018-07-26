/*
  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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

#include <python.hpp>
#include "Interaction.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(Interaction::theLogger, "Interaction");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////

    void
    Interaction::registerPython() {
      using namespace espressopp::python;

      real (Interaction::*pyComputeEnergyAAraw)() = &Interaction::computeEnergyAA;
      real (Interaction::*pyComputeEnergyAAtype)(int) = &Interaction::computeEnergyAA;
      real (Interaction::*pyComputeEnergyCGraw)() = &Interaction::computeEnergyCG;
      real (Interaction::*pyComputeEnergyCGtype)(int) = &Interaction::computeEnergyCG;

      class_< Interaction, boost::noncopyable >("interaction_Interaction", no_init)
        .def("computeEnergy", &Interaction::computeEnergy)
        .def("computeEnergyDeriv", &Interaction::computeEnergyDeriv)
        .def("computeEnergyAA", pyComputeEnergyAAraw)
        .def("computeEnergyAA", pyComputeEnergyAAtype)
        .def("computeEnergyCG", pyComputeEnergyCGraw)
        .def("computeEnergyCG", pyComputeEnergyCGtype)
        .def("computeVirial", &Interaction::computeVirial)
        .def("bondType", &Interaction::bondType)
      ;
    }
  }
}
