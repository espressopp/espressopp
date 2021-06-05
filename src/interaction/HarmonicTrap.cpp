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
#include "HarmonicTrap.hpp"
#include "Real3D.hpp"

namespace espressopp
{
namespace interaction
{
typedef class SingleParticleInteractionTemplate<HarmonicTrap> SingleParticleHarmonicTrap;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void HarmonicTrap::registerPython()
{
    using namespace espressopp::python;

    class_<HarmonicTrap, bases<SingleParticlePotential> >("interaction_HarmonicTrap", init<>())
        .add_property("k", &HarmonicTrap::getK, &HarmonicTrap::setK)
        .add_property("center", &HarmonicTrap::getCenter, &HarmonicTrap::setCenter);

    class_<SingleParticleHarmonicTrap, bases<Interaction> >(
        "interaction_SingleParticleHarmonicTrap",
        init<std::shared_ptr<System>, std::shared_ptr<HarmonicTrap> >())
        .def("setPotential", &SingleParticleHarmonicTrap::setPotential)
        .def("getPotential", &SingleParticleHarmonicTrap::getPotential);
}
}  // namespace interaction
}  // namespace espressopp
