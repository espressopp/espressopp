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
#include "SquareWell.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate < SquareWell > VerletListSquareWell;
    typedef class FixedPairListInteractionTemplate< SquareWell > FixedPairListSquareWell;
    typedef class FixedPairListTypesInteractionTemplate< SquareWell > FixedPairTypesSquareWell;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void SquareWell::registerPython() {
      using namespace espressopp::python;

      class_ <SquareWell, bases <Potential> >
    	("interaction_SquareWell", init< real, real, real >())
        .def(init< real, real, real, real >())
        .add_property("epsilon", &SquareWell::getEpsilon, &SquareWell::setEpsilon)
        .add_property("sigma", &SquareWell::getSigma, &SquareWell::setSigma)
        .add_property("width", &SquareWell::getLambda, &SquareWell::setLambda)
        .add_property("a", &SquareWell::getA, &SquareWell::setA)
        .def_pickle(SquareWell_pickle())
        ;

      class_ <VerletListSquareWell, bases <Interaction> >
        ("interaction_VerletListSquareWell", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListSquareWell::getVerletList)
        .def("setPotential", &VerletListSquareWell::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListSquareWell::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_ <FixedPairListSquareWell, bases <Interaction> >
        ("interaction_FixedPairListSquareWell",
         init< shared_ptr<System>, shared_ptr<FixedPairList>,  shared_ptr<SquareWell> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<SquareWell> >())
        .def("setPotential", &FixedPairListSquareWell::setPotential)
        .def("setFixedPairList", &FixedPairListSquareWell::setFixedPairList)
        .def("getFixedPairList", &FixedPairListSquareWell::getFixedPairList)
        ;
    }
  }
}
