/*
  Copyright (C) 2014
      Pierre de Buyl
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
#include "SingleParticlePotential.hpp"
#include "bc/BC.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(SingleParticlePotential::theLogger, "SingleParticlePotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void SingleParticlePotential::registerPython() {
        using namespace espressopp::python;

        real (SingleParticlePotential::*computeEnergy)(const Particle& p, const bc::BC& bc) const =
          &SingleParticlePotential::computeEnergy;

        Real3D (SingleParticlePotential::*computeForce)(const Particle& p, const bc::BC& bc) const =
          &SingleParticlePotential::computeForce;

        class_< SingleParticlePotential, boost::noncopyable >
            ("interaction_SingleParticlePotential", no_init)
            .def("computeEnergy", pure_virtual(computeEnergy))
            .def("computeForce", pure_virtual(computeForce))
        ;
    }
  }
}
