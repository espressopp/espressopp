/*
  Copyright (C) 2019-2021
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
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
#include "LennardJones.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListLennardJones.hpp"

namespace espressopp
{
namespace vec
{
namespace interaction
{
LOG4ESPP_LOGGER(LennardJones::theLogger, "LennardJones");

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void LennardJones::registerPython()
{
    using namespace espressopp::python;

    class_<LennardJones, bases<Potential> >("vec_interaction_LennardJones",
                                            init<real, real, real>())
        .def(init<real, real, real, real>())
        .add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
        .add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
        .def_pickle(LennardJones_pickle());

    class_<VerletListLennardJonesBase, bases<Interaction> >(
        "vec_interaction_VerletListLennardJonesBase", init<std::shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesBase::getVerletList)
        .def("setPotential", &VerletListLennardJonesBase::setPotential)
        .def("getPotential", &VerletListLennardJonesBase::getPotentialPtr);

    class_<VerletListLennardJones, bases<Interaction> >("vec_interaction_VerletListLennardJones",
                                                        init<std::shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJones::getVerletList)
        .def("setPotential", &VerletListLennardJones::setPotential)
        .def("getPotential", &VerletListLennardJones::getPotentialPtr);
}
}  // namespace interaction
}  // namespace vec
}  // namespace espressopp
