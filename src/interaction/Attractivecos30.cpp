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
#include "Attractivecos30.hpp"
#include "Tabulated.hpp"
#include "Harmonic.hpp"
#include "ReactionFieldGeneralized.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListAdressATInteractionTemplate.hpp"
#include "VerletListAdressCGInteractionTemplate.hpp"
#include "VerletListAdressATATInteractionTemplate.hpp"
#include "VerletListAdressATATCGInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "VerletListHadressATInteractionTemplate.hpp"
#include "VerletListHadressCGInteractionTemplate.hpp"
#include "VerletListHadressATATInteractionTemplate.hpp"
#include "VerletListHadressATATCGInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <Attractivecos30>
        VerletListAttractivecos30;
    typedef class VerletListAdressInteractionTemplate <Attractivecos30, Tabulated>
        VerletListAdressAttractivecos30;
    typedef class VerletListAdressATInteractionTemplate <Attractivecos30>
        VerletListAdressATAttractivecos30;
    typedef class VerletListAdressATATInteractionTemplate <Attractivecos30, ReactionFieldGeneralized>
        VerletListAdressATLenJonesReacFieldGen;
    typedef class VerletListAdressATATCGInteractionTemplate <Attractivecos30, ReactionFieldGeneralized, Tabulated>
        VerletListAdressATLJReacFieldGenTab;
    typedef class VerletListAdressATATCGInteractionTemplate <Attractivecos30, ReactionFieldGeneralized, Harmonic>
        VerletListAdressATLJReacFieldGenHarmonic;
    typedef class VerletListAdressCGInteractionTemplate <Attractivecos30>
        VerletListAdressCGAttractivecos30;
    typedef class VerletListAdressInteractionTemplate <Attractivecos30, Attractivecos30>
        VerletListAdressAttractivecos302;
    typedef class VerletListAdressInteractionTemplate <Attractivecos30, Harmonic>
        VerletListAdressAttractivecos30Harmonic;
    typedef class VerletListHadressInteractionTemplate <Attractivecos30, Tabulated>
        VerletListHadressAttractivecos30;
    typedef class VerletListHadressATInteractionTemplate <Attractivecos30>
        VerletListHadressATAttractivecos30;
    typedef class VerletListHadressATATInteractionTemplate <Attractivecos30, ReactionFieldGeneralized>
        VerletListHadressATLenJonesReacFieldGen;
    typedef class VerletListHadressATATCGInteractionTemplate <Attractivecos30, ReactionFieldGeneralized, Tabulated>
        VerletListHadressATLJReacFieldGenTab;
    typedef class VerletListHadressATATCGInteractionTemplate <Attractivecos30, ReactionFieldGeneralized, Harmonic>
        VerletListHadressATLJReacFieldGenHarmonic;
    typedef class VerletListHadressCGInteractionTemplate <Attractivecos30>
        VerletListHadressCGAttractivecos30;
    typedef class VerletListHadressInteractionTemplate <Attractivecos30, Attractivecos30>
        VerletListHadressAttractivecos302;
    typedef class VerletListHadressInteractionTemplate <Attractivecos30, Harmonic>
        VerletListHadressAttractivecos30Harmonic;
    typedef class CellListAllPairsInteractionTemplate <Attractivecos30>
        CellListAttractivecos30;
    typedef class FixedPairListInteractionTemplate <Attractivecos30>
        FixedPairListAttractivecos30;
    typedef class FixedPairListTypesInteractionTemplate <Attractivecos30>
        FixedPairListTypesAttractivecos30;
    LOG4ESPP_LOGGER(Attractivecos30::theLogger, "Attractivecos30");
    // LOG4ESPP_LOGGER(VerletListAttractivecos30::theLogger, "VerletListAttractivecos30");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Attractivecos30::registerPython() {
      using namespace espressopp::python;

      class_< Attractivecos30, bases< Potential > >
        ("interaction_Attractivecos30", init< real, real, real >())
            .def(init< real, real, real, real >())
        .add_property("mu", &Attractivecos30::getMu, &Attractivecos30::setMu)
        .add_property("alpha", &Attractivecos30::getAlpha, &Attractivecos30::setAlpha)
        .def_pickle(Attractivecos30_pickle())

      ;

      class_< VerletListAttractivecos30, bases< Interaction > >
        ("interaction_VerletListAttractivecos30", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListAttractivecos30::getVerletList)
        .def("setPotential", &VerletListAttractivecos30::setPotential)
        .def("getPotential", &VerletListAttractivecos30::getPotentialPtr)
      ;

      class_< VerletListAdressATAttractivecos30, bases< Interaction > >
        ("interaction_VerletListAdressATAttractivecos30",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressATAttractivecos30::getVerletList)
        .def("setPotential", &VerletListAdressATAttractivecos30::setPotential)
        .def("getPotential", &VerletListAdressATAttractivecos30::getPotentialPtr)
      ;

      class_< VerletListAdressATLenJonesReacFieldGen, bases< Interaction > >
        ("interaction_VerletListAdressATLenJonesReacFieldGen",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotential1", &VerletListAdressATLenJonesReacFieldGen::setPotential1)
        .def("setPotential2", &VerletListAdressATLenJonesReacFieldGen::setPotential2);
      ;

      class_< VerletListAdressATLJReacFieldGenTab, bases< Interaction > >
        ("interaction_VerletListAdressATLJReacFieldGenTab",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListAdressATLJReacFieldGenTab::setPotentialAT1)
        .def("setPotentialAT2", &VerletListAdressATLJReacFieldGenTab::setPotentialAT2)
        .def("setPotentialCG", &VerletListAdressATLJReacFieldGenTab::setPotentialCG);
      ;

      class_< VerletListAdressATLJReacFieldGenHarmonic, bases< Interaction > >
        ("interaction_VerletListAdressATLJReacFieldGenHarmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListAdressATLJReacFieldGenHarmonic::setPotentialAT1)
        .def("setPotentialAT2", &VerletListAdressATLJReacFieldGenHarmonic::setPotentialAT2)
        .def("setPotentialCG", &VerletListAdressATLJReacFieldGenHarmonic::setPotentialCG);
      ;

      class_< VerletListAdressCGAttractivecos30, bases< Interaction > >
        ("interaction_VerletListAdressCGAttractivecos30",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressCGAttractivecos30::getVerletList)
        .def("setPotential", &VerletListAdressCGAttractivecos30::setPotential)
        .def("getPotential", &VerletListAdressCGAttractivecos30::getPotentialPtr)
      ;

      class_< VerletListAdressAttractivecos30, bases< Interaction > >
        ("interaction_VerletListAdressAttractivecos30",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressAttractivecos30::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressAttractivecos30::setPotentialCG);
      ;

      class_< VerletListAdressAttractivecos302, bases< Interaction > >
        ("interaction_VerletListAdressAttractivecos302",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressAttractivecos302::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressAttractivecos302::setPotentialCG);
      ;

      class_< VerletListAdressAttractivecos30Harmonic, bases< Interaction > >
        ("interaction_VerletListAdressAttractivecos30Harmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressAttractivecos30Harmonic::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressAttractivecos30Harmonic::setPotentialCG);
      ;

      class_< VerletListHadressATAttractivecos30, bases< Interaction > >
        ("interaction_VerletListHadressATAttractivecos30",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressATAttractivecos30::getVerletList)
        .def("setPotential", &VerletListHadressATAttractivecos30::setPotential)
        .def("getPotential", &VerletListHadressATAttractivecos30::getPotentialPtr)
      ;

      class_< VerletListHadressATLenJonesReacFieldGen, bases< Interaction > >
        ("interaction_VerletListHadressATLenJonesReacFieldGen",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotential1", &VerletListHadressATLenJonesReacFieldGen::setPotential1)
        .def("setPotential2", &VerletListHadressATLenJonesReacFieldGen::setPotential2);
      ;

      class_< VerletListHadressATLJReacFieldGenTab, bases< Interaction > >
        ("interaction_VerletListHadressATLJReacFieldGenTab",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListHadressATLJReacFieldGenTab::setPotentialAT1)
        .def("setPotentialAT2", &VerletListHadressATLJReacFieldGenTab::setPotentialAT2)
        .def("setPotentialCG", &VerletListHadressATLJReacFieldGenTab::setPotentialCG);
      ;

      class_< VerletListHadressATLJReacFieldGenHarmonic, bases< Interaction > >
        ("interaction_VerletListHadressATLJReacFieldGenHarmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListHadressATLJReacFieldGenHarmonic::setPotentialAT1)
        .def("setPotentialAT2", &VerletListHadressATLJReacFieldGenHarmonic::setPotentialAT2)
        .def("setPotentialCG", &VerletListHadressATLJReacFieldGenHarmonic::setPotentialCG);
      ;

      class_< VerletListHadressCGAttractivecos30, bases< Interaction > >
        ("interaction_VerletListHadressCGAttractivecos30",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressCGAttractivecos30::getVerletList)
        .def("setPotential", &VerletListHadressCGAttractivecos30::setPotential)
        .def("getPotential", &VerletListHadressCGAttractivecos30::getPotentialPtr)
      ;

      class_< VerletListHadressAttractivecos30, bases< Interaction > >
        ("interaction_VerletListHadressAttractivecos30",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressAttractivecos30::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressAttractivecos30::setPotentialCG);
      ;

      class_< VerletListHadressAttractivecos302, bases< Interaction > >
        ("interaction_VerletListHadressAttractivecos302",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressAttractivecos302::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressAttractivecos302::setPotentialCG);
      ;

      class_< VerletListHadressAttractivecos30Harmonic, bases< Interaction > >
        ("interaction_VerletListHadressAttractivecos30Harmonic",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressAttractivecos30Harmonic::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressAttractivecos30Harmonic::setPotentialCG);
      ;

      class_< CellListAttractivecos30, bases< Interaction > >
        ("interaction_CellListAttractivecos30", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListAttractivecos30::setPotential);
    ;

      class_< FixedPairListAttractivecos30, bases< Interaction > >
        ("interaction_FixedPairListAttractivecos30",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Attractivecos30> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Attractivecos30> >())
          .def("setPotential", &FixedPairListAttractivecos30::setPotential)
          .def("getPotential", &FixedPairListAttractivecos30::getPotential)
          .def("setFixedPairList", &FixedPairListAttractivecos30::setFixedPairList)
          .def("getFixedPairList", &FixedPairListAttractivecos30::getFixedPairList)
      ;

      class_< FixedPairListTypesAttractivecos30, bases< Interaction > >
        ("interaction_FixedPairListTypesAttractivecos30", init< shared_ptr<System>, shared_ptr<FixedPairList> >())
         .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
         .def("setFixedPairList", &FixedPairListTypesAttractivecos30::setFixedPairList)
         .def("getFixedPairList", &FixedPairListTypesAttractivecos30::getFixedPairList)
         .def("setPotential", &FixedPairListTypesAttractivecos30::setPotential)
         .def("getPotential", &FixedPairListTypesAttractivecos30::getPotentialPtr)
      ;
    }

  }
}
