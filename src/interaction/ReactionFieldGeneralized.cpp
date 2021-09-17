/*
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
#include "ReactionFieldGeneralized.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListAdressATInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "VerletListHadressATInteractionTemplate.hpp"
#include "VerletListHybridInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"

namespace espressopp
{
namespace interaction
{
typedef class VerletListInteractionTemplate<ReactionFieldGeneralized>
    VerletListReactionFieldGeneralized;

typedef class VerletListHybridInteractionTemplate<ReactionFieldGeneralized>
    VerletListHybridReactionFieldGeneralized;

typedef class VerletListAdressInteractionTemplate<ReactionFieldGeneralized, Tabulated>
    VerletListAdressReactionFieldGeneralized;

typedef class VerletListAdressATInteractionTemplate<ReactionFieldGeneralized>
    VerletListAdressATReactionFieldGeneralized;

typedef class VerletListHadressInteractionTemplate<ReactionFieldGeneralized, Tabulated>
    VerletListHadressReactionFieldGeneralized;

typedef class VerletListHadressATInteractionTemplate<ReactionFieldGeneralized>
    VerletListHadressATReactionFieldGeneralized;

typedef class CellListAllPairsInteractionTemplate<ReactionFieldGeneralized>
    CellListReactionFieldGeneralized;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void ReactionFieldGeneralized::registerPython()
{
    using namespace espressopp::python;

    class_<ReactionFieldGeneralized, bases<Potential> >("interaction_ReactionFieldGeneralized",
                                                        init<real, real, real, real, real>())
        .def(init<real, real, real, real, real, real>())
        .def_pickle(ReactionFieldGeneralized_pickle())
        .add_property("prefactor", &ReactionFieldGeneralized::getPrefactor,
                      &ReactionFieldGeneralized::setPrefactor);

    class_<VerletListReactionFieldGeneralized, bases<Interaction> >(
        "interaction_VerletListReactionFieldGeneralized", init<std::shared_ptr<VerletList> >())
        .def("setPotential", &VerletListReactionFieldGeneralized::setPotential)
        .def("getPotential", &VerletListReactionFieldGeneralized::getPotentialPtr);

    class_<VerletListAdressATReactionFieldGeneralized, bases<Interaction> >(
        "interaction_VerletListAdressATReactionFieldGeneralized",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotential", &VerletListAdressATReactionFieldGeneralized::setPotential)
        .def("getPotential", &VerletListAdressATReactionFieldGeneralized::getPotentialPtr);

    class_<VerletListHybridReactionFieldGeneralized, bases<Interaction> >(
        "interaction_VerletListHybridReactionFieldGeneralized",
        init<shared_ptr<VerletList>, bool>())
        .def("setPotential", &VerletListHybridReactionFieldGeneralized::setPotential)
        .def("getPotential", &VerletListHybridReactionFieldGeneralized::getPotentialPtr)
        .add_property("scale_factor", &VerletListHybridReactionFieldGeneralized::scaleFactor,
                      &VerletListHybridReactionFieldGeneralized::setScaleFactor)
        .add_property("max_force", &VerletListHybridReactionFieldGeneralized::maxForce,
                      &VerletListHybridReactionFieldGeneralized::setMaxForce);

    class_<VerletListAdressReactionFieldGeneralized, bases<Interaction> >(
        "interaction_VerletListAdressReactionFieldGeneralized",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressReactionFieldGeneralized::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressReactionFieldGeneralized::setPotentialCG);
    ;

    class_<VerletListHadressATReactionFieldGeneralized, bases<Interaction> >(
        "interaction_VerletListHadressATReactionFieldGeneralized",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotential", &VerletListHadressATReactionFieldGeneralized::setPotential)
        .def("getPotential", &VerletListHadressATReactionFieldGeneralized::getPotentialPtr);

    class_<VerletListHadressReactionFieldGeneralized, bases<Interaction> >(
        "interaction_VerletListHadressReactionFieldGeneralized",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressReactionFieldGeneralized::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressReactionFieldGeneralized::setPotentialCG);
    ;

    class_<CellListReactionFieldGeneralized, bases<Interaction> >(
        "interaction_CellListReactionFieldGeneralized", init<std::shared_ptr<storage::Storage> >())
        .def("setPotential", &CellListReactionFieldGeneralized::setPotential);
    ;
}
}  // namespace interaction
}  // namespace espressopp
