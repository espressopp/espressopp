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
#include "Harmonic.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressATInteractionTemplate.hpp"
#include "VerletListAdressCGInteractionTemplate.hpp"
#include "VerletListHadressATInteractionTemplate.hpp"
#include "VerletListHadressCGInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Harmonic::registerPython() {
      using namespace espressopp::python;

      class_< Harmonic, bases< Potential > >
    	("interaction_Harmonic", init< real, real, real >())
	.def(init< real, real, real, real >())
	.add_property("K", &Harmonic::getK, &Harmonic::setK)
	.add_property("r0", &Harmonic::getR0, &Harmonic::setR0)
    	;

      typedef class FixedPairListInteractionTemplate< Harmonic >
        FixedPairListHarmonic;
      typedef class FixedPairListTypesInteractionTemplate< Harmonic >
        FixedPairListTypesHarmonic;
      typedef class VerletListHadressATInteractionTemplate< Harmonic >
        VerletListHadressATHarmonic;
      typedef class VerletListHadressCGInteractionTemplate< Harmonic >
        VerletListHadressCGHarmonic;
      typedef class VerletListAdressATInteractionTemplate< Harmonic >
        VerletListAdressATHarmonic;
      typedef class VerletListAdressCGInteractionTemplate< Harmonic >
        VerletListAdressCGHarmonic;
      typedef class VerletListInteractionTemplate< Harmonic >
        VerletListHarmonic;

      class_< FixedPairListHarmonic, bases< Interaction > >
        ("interaction_FixedPairListHarmonic",
           init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Harmonic> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Harmonic> >())
        .def("setPotential", &FixedPairListHarmonic::setPotential)
        .def("getPotential", &FixedPairListHarmonic::getPotential)
        .def("setFixedPairList", &FixedPairListHarmonic::setFixedPairList)
        .def("getFixedPairList", &FixedPairListHarmonic::getFixedPairList);
      ;
      class_< FixedPairListTypesHarmonic, bases< Interaction > >
        ("interaction_FixedPairListTypesHarmonic",
           init< shared_ptr<System>, shared_ptr<FixedPairList> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
        .def("setPotential", &FixedPairListTypesHarmonic::setPotential)
        .def("getPotential", &FixedPairListTypesHarmonic::getPotentialPtr)
        .def("setFixedPairList", &FixedPairListTypesHarmonic::setFixedPairList)
        .def("getFixedPairList", &FixedPairListTypesHarmonic::getFixedPairList);
     ;

      class_< VerletListHarmonic, bases< Interaction > >
        ("interaction_VerletListHarmonic",
           init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListHarmonic::setPotential)
        .def("getPotential", &VerletListHarmonic::getPotentialPtr)
      ;

            class_< VerletListAdressATHarmonic, bases< Interaction > >
        ("interaction_VerletListAdressATHarmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressATHarmonic::getVerletList)
        .def("setPotential", &VerletListAdressATHarmonic::setPotential)
        .def("getPotential", &VerletListAdressATHarmonic::getPotentialPtr)
      ;

      class_< VerletListAdressCGHarmonic, bases< Interaction > >
        ("interaction_VerletListAdressCGHarmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressCGHarmonic::getVerletList)
        .def("setPotential", &VerletListAdressCGHarmonic::setPotential)
        .def("getPotential", &VerletListAdressCGHarmonic::getPotentialPtr)
      ;

            class_< VerletListHadressATHarmonic, bases< Interaction > >
        ("interaction_VerletListHadressATHarmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressATHarmonic::getVerletList)
        .def("setPotential", &VerletListHadressATHarmonic::setPotential)
        .def("getPotential", &VerletListHadressATHarmonic::getPotentialPtr)
      ;

      class_< VerletListHadressCGHarmonic, bases< Interaction > >
        ("interaction_VerletListHadressCGHarmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressCGHarmonic::getVerletList)
        .def("setPotential", &VerletListHadressCGHarmonic::setPotential)
        .def("getPotential", &VerletListHadressCGHarmonic::getPotentialPtr)
      ;
    }

  }
}
