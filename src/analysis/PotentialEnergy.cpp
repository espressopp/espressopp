/*
  Copyright (C) 2015
    Jakub Krajniak (jkrajniak at gmail.com)

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
#include "PotentialEnergy.hpp"
#include "interaction/Interaction.hpp"

using namespace espressopp;  //NOLINT

namespace espressopp {
namespace analysis {

real PotentialEnergy::compute_real() const {
  if (compute_global_)
    return interaction_->computeEnergy();
  else if (compute_at_)
    return interaction_->computeEnergyAA();
  else
    return interaction_->computeEnergyCG();
}


void PotentialEnergy::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<PotentialEnergy, bases<Observable> >
    ("analysis_PotentialEnergy",
        init< shared_ptr<System>, shared_ptr<interaction::Interaction> >())
    .def(init<shared_ptr<System>, shared_ptr<interaction::Interaction>, bool>())
    .add_property("value", &PotentialEnergy::compute_real);
}
}  // end namespace analysis
}  // end namespace espressopp
