/*
  Copyright (C) 2017,2018
  Max Planck Institute for Polymer Research

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
#include "SmoothSquareWell.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate < SmoothSquareWell > VerletListSmoothSquareWell;
    typedef class FixedPairListInteractionTemplate< SmoothSquareWell > FixedPairListSmoothSquareWell;
    typedef class FixedPairListTypesInteractionTemplate< SmoothSquareWell > FixedPairListTypesSmoothSquareWell;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void SmoothSquareWell::registerPython() {
      using namespace espressopp::python;

      class_ <SmoothSquareWell, bases <Potential> >
        ("interaction_SmoothSquareWell", init< real, real, real >())
        .def(init< real, real, real, real >())
        .add_property("epsilon", &SmoothSquareWell::getEpsilon, &SmoothSquareWell::setEpsilon)
        .add_property("sigma", &SmoothSquareWell::getSigma, &SmoothSquareWell::setSigma)
        .add_property("Lambda", &SmoothSquareWell::getLambda, &SmoothSquareWell::setLambda)
        .add_property("a", &SmoothSquareWell::getA, &SmoothSquareWell::setA)
        .def_pickle(SmoothSquareWell_pickle())
        ;

      class_ <VerletListSmoothSquareWell, bases <Interaction> >
        ("interaction_VerletListSmoothSquareWell", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListSmoothSquareWell::getVerletList)
        .def("setPotential", &VerletListSmoothSquareWell::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListSmoothSquareWell::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_ <FixedPairListSmoothSquareWell, bases <Interaction> >
        ("interaction_FixedPairListSmoothSquareWell",
         init< shared_ptr<System>, shared_ptr<FixedPairList>,  shared_ptr<SmoothSquareWell> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<SmoothSquareWell> >())
        .def("setPotential", &FixedPairListSmoothSquareWell::setPotential)
        .def("getPotential", &FixedPairListSmoothSquareWell::getPotential)
        .def("setFixedPairList", &FixedPairListSmoothSquareWell::setFixedPairList)
        .def("getFixedPairList", &FixedPairListSmoothSquareWell::getFixedPairList)
        ;

      class_ <FixedPairListTypesSmoothSquareWell, bases <Interaction> >
        ("interaction_FixedPairListTypesSmoothSquareWell",
         init< shared_ptr<System>, shared_ptr<FixedPairList> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
        .def("setPotential", &FixedPairListTypesSmoothSquareWell::setPotential)
        .def("getPotential", &FixedPairListTypesSmoothSquareWell::getPotentialPtr)
        .def("setFixedPairList", &FixedPairListTypesSmoothSquareWell::setFixedPairList)
        .def("getFixedPairList", &FixedPairListTypesSmoothSquareWell::getFixedPairList)
        ;
    }
  }
}
