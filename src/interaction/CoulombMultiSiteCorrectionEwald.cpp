/*
  Copyright (C) 2022
      Data Center, Johannes Gutenberg University Mainz

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
#include "CoulombMultiSiteCorrectionEwald.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

// For a Coulombic FixedPairList interaction, it's necessary to use
// FixedPairListTypesInteractionTemplate.hpp instead of FixedPairListInteractionTemplate.hpp so that
// we can use _computeForce and _computeEnergy which take both particles and distance vector as
// arguments because the Coulomb interaction needs access to both the charges (via the particles)
// and the minimum image distance (via the boundary conditions in the interaction template)

namespace espressopp
{
namespace interaction
{
typedef class VerletListInteractionTemplate<CoulombMultiSiteCorrectionEwald>
    VerletListCoulombMultiSiteCorrectionEwald;
typedef class FixedPairListTypesInteractionTemplate<CoulombMultiSiteCorrectionEwald>
    FixedPairListTypesCoulombMultiSiteCorrectionEwald;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void CoulombMultiSiteCorrectionEwald::registerPython()
{
    using namespace espressopp::python;

    class_<CoulombMultiSiteCorrectionEwald, bases<Potential> >(
        "interaction_CoulombMultiSiteCorrectionEwald", init<>())
        .def(init<real, real, real>())
        .add_property("alpha", &CoulombMultiSiteCorrectionEwald::getAlpha,
                      &CoulombMultiSiteCorrectionEwald::setAlpha)
        .add_property("prefactor", &CoulombMultiSiteCorrectionEwald::getPrefactor,
                      &CoulombMultiSiteCorrectionEwald::setPrefactor);

    class_<VerletListCoulombMultiSiteCorrectionEwald, bases<Interaction> >(
        "interaction_VerletListCoulombMultiSiteCorrectionEwald", init<shared_ptr<VerletList> >())
        .def("setPotential", &VerletListCoulombMultiSiteCorrectionEwald::setPotential)
        .def("getPotential", &VerletListCoulombMultiSiteCorrectionEwald::getPotentialPtr);

    class_<FixedPairListTypesCoulombMultiSiteCorrectionEwald, bases<Interaction> >(
        "interaction_FixedPairListTypesCoulombMultiSiteCorrectionEwald",
        init<shared_ptr<System>, shared_ptr<FixedPairList> >())
        .def(init<shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
        .def("setPotential", &FixedPairListTypesCoulombMultiSiteCorrectionEwald::setPotential);
}

}  // namespace interaction
}  // namespace espressopp
