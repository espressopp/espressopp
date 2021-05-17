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
#include "LennardJonesAutoBonds.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp
{
namespace interaction
{
typedef class VerletListInteractionTemplate<LennardJonesAutoBonds> VerletListLennardJonesAutoBonds;
typedef class VerletListAdressInteractionTemplate<LennardJonesAutoBonds, Tabulated>
    VerletListAdressLennardJonesAutoBonds;
typedef class VerletListHadressInteractionTemplate<LennardJonesAutoBonds, Tabulated>
    VerletListHadressLennardJonesAutoBonds;
typedef class CellListAllPairsInteractionTemplate<LennardJonesAutoBonds>
    CellListLennardJonesAutoBonds;
typedef class FixedPairListInteractionTemplate<LennardJonesAutoBonds>
    FixedPairListLennardJonesAutoBonds;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void LennardJonesAutoBonds::registerPython()
{
    using namespace espressopp::python;

    class_<LennardJonesAutoBonds, bases<Potential> >(
        "interaction_LennardJonesAutoBonds",
        init<real, real, real, std::shared_ptr<FixedPairList>, int>())
        .def(init<real, real, real, real, std::shared_ptr<FixedPairList>, int>())
        .add_property("sigma", &LennardJonesAutoBonds::getSigma, &LennardJonesAutoBonds::setSigma)
        .add_property("epsilon", &LennardJonesAutoBonds::getEpsilon,
                      &LennardJonesAutoBonds::setEpsilon)
        .add_property("max_crosslinks", &LennardJonesAutoBonds::getMaxCrosslinks,
                      &LennardJonesAutoBonds::setMaxCrosslinks);

    class_<VerletListLennardJonesAutoBonds, bases<Interaction> >(
        "interaction_VerletListLennardJonesAutoBonds", init<std::shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesAutoBonds::getVerletList)
        .def("setPotential", &VerletListLennardJonesAutoBonds::setPotential,
             return_value_policy<reference_existing_object>())
        .def("getPotential", &VerletListLennardJonesAutoBonds::getPotential,
             return_value_policy<reference_existing_object>());

    class_<VerletListAdressLennardJonesAutoBonds, bases<Interaction> >(
        "interaction_VerletListAdressLennardJonesAutoBonds",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesAutoBonds::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesAutoBonds::setPotentialCG);
    ;

    class_<VerletListHadressLennardJonesAutoBonds, bases<Interaction> >(
        "interaction_VerletListHadressLennardJonesAutoBonds",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJonesAutoBonds::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJonesAutoBonds::setPotentialCG);
    ;

    class_<CellListLennardJonesAutoBonds, bases<Interaction> >(
        "interaction_CellListLennardJonesAutoBonds", init<std::shared_ptr<storage::Storage> >())
        .def("setPotential", &CellListLennardJonesAutoBonds::setPotential);
    ;

    class_<FixedPairListLennardJonesAutoBonds, bases<Interaction> >(
        "interaction_FixedPairListLennardJonesAutoBonds",
        init<std::shared_ptr<System>, std::shared_ptr<FixedPairList>,
             std::shared_ptr<LennardJonesAutoBonds> >())
        .def(init<std::shared_ptr<System>, std::shared_ptr<FixedPairListAdress>,
                  std::shared_ptr<LennardJonesAutoBonds> >())
        .def("setPotential", &FixedPairListLennardJonesAutoBonds::setPotential);
    ;
}

}  // namespace interaction
}  // namespace espressopp
