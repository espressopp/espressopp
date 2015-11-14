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
#include "LennardJones93Wall.hpp"
#include "Real3D.hpp"

namespace espressopp {
  namespace interaction {

    typedef class SingleParticleInteractionTemplate <LennardJones93Wall>
    SingleParticleLennardJones93Wall;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    LennardJones93Wall::registerPython() {
      using namespace espressopp::python;

      class_< LennardJones93Wall, bases< SingleParticlePotential > >
        ("interaction_LennardJones93Wall", init<>())
        .def("setParams", &LennardJones93Wall::setParams)
        .def("getParams", &LennardJones93Wall::getParams)
        ;

      class_< SingleParticleLennardJones93Wall, bases< Interaction > >
        ("interaction_SingleParticleLennardJones93Wall", init< shared_ptr<System>, shared_ptr<LennardJones93Wall> >())
        .def("setPotential", &SingleParticleLennardJones93Wall::setPotential)
        .def("getPotential", &SingleParticleLennardJones93Wall::getPotential)
       ;
    }
  }
}
