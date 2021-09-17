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
#include "LennardJones.hpp"
#include "Tabulated.hpp"
#include "Harmonic.hpp"
#include "ReactionFieldGeneralized.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListHybridInteractionTemplate.hpp"
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
#include "FixedPairListAdressInteractionTemplate.hpp"

namespace espressopp
{
namespace interaction
{
typedef class VerletListInteractionTemplate<LennardJones> VerletListLennardJones;
typedef class VerletListHybridInteractionTemplate<LennardJones> VerletListHybridLennardJones;
typedef class VerletListAdressInteractionTemplate<LennardJones, Tabulated>
    VerletListAdressLennardJones;
typedef class VerletListAdressATInteractionTemplate<LennardJones> VerletListAdressATLennardJones;
typedef class VerletListAdressATATInteractionTemplate<LennardJones, ReactionFieldGeneralized>
    VerletListAdressATLenJonesReacFieldGen;
typedef class VerletListAdressATATCGInteractionTemplate<LennardJones,
                                                        ReactionFieldGeneralized,
                                                        Tabulated>
    VerletListAdressATLJReacFieldGenTab;
typedef class VerletListAdressATATCGInteractionTemplate<LennardJones,
                                                        ReactionFieldGeneralized,
                                                        Harmonic>
    VerletListAdressATLJReacFieldGenHarmonic;
typedef class VerletListAdressCGInteractionTemplate<LennardJones> VerletListAdressCGLennardJones;
typedef class VerletListAdressInteractionTemplate<LennardJones, LennardJones>
    VerletListAdressLennardJones2;
typedef class VerletListAdressInteractionTemplate<LennardJones, Harmonic>
    VerletListAdressLennardJonesHarmonic;
typedef class VerletListHadressInteractionTemplate<LennardJones, Tabulated>
    VerletListHadressLennardJones;
typedef class VerletListHadressATInteractionTemplate<LennardJones> VerletListHadressATLennardJones;
typedef class VerletListHadressATATInteractionTemplate<LennardJones, ReactionFieldGeneralized>
    VerletListHadressATLenJonesReacFieldGen;
typedef class VerletListHadressATATCGInteractionTemplate<LennardJones,
                                                         ReactionFieldGeneralized,
                                                         Tabulated>
    VerletListHadressATLJReacFieldGenTab;
typedef class VerletListHadressATATCGInteractionTemplate<LennardJones,
                                                         ReactionFieldGeneralized,
                                                         Harmonic>
    VerletListHadressATLJReacFieldGenHarmonic;
typedef class VerletListHadressCGInteractionTemplate<LennardJones> VerletListHadressCGLennardJones;
typedef class VerletListHadressInteractionTemplate<LennardJones, LennardJones>
    VerletListHadressLennardJones2;
typedef class VerletListHadressInteractionTemplate<LennardJones, Harmonic>
    VerletListHadressLennardJonesHarmonic;
typedef class CellListAllPairsInteractionTemplate<LennardJones> CellListLennardJones;
typedef class FixedPairListInteractionTemplate<LennardJones> FixedPairListLennardJones;
typedef class FixedPairListTypesInteractionTemplate<LennardJones> FixedPairListTypesLennardJones;
typedef class FixedPairListAdressInteractionTemplate<LennardJones> FixedPairListAdressLennardJones;
LOG4ESPP_LOGGER(LennardJones::theLogger, "LennardJones");

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void LennardJones::registerPython()
{
    using namespace espressopp::python;

    class_<LennardJones, bases<Potential> >("interaction_LennardJones", init<real, real, real>())
        .def(init<real, real, real, real>())
        .add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
        .add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
        .def_pickle(LennardJones_pickle())

        ;

    class_<VerletListLennardJones, bases<Interaction> >("interaction_VerletListLennardJones",
                                                        init<std::shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJones::getVerletList)
        .def("setPotential", &VerletListLennardJones::setPotential)
        .def("getPotential", &VerletListLennardJones::getPotentialPtr);

    class_<VerletListAdressATLennardJones, bases<Interaction> >(
        "interaction_VerletListAdressATLennardJones",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressATLennardJones::getVerletList)
        .def("setPotential", &VerletListAdressATLennardJones::setPotential)
        .def("getPotential", &VerletListAdressATLennardJones::getPotentialPtr);

    class_<VerletListAdressATLenJonesReacFieldGen, bases<Interaction> >(
        "interaction_VerletListAdressATLenJonesReacFieldGen",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotential1", &VerletListAdressATLenJonesReacFieldGen::setPotential1)
        .def("setPotential2", &VerletListAdressATLenJonesReacFieldGen::setPotential2);
    ;

    class_<VerletListAdressATLJReacFieldGenTab, bases<Interaction> >(
        "interaction_VerletListAdressATLJReacFieldGenTab",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListAdressATLJReacFieldGenTab::setPotentialAT1)
        .def("setPotentialAT2", &VerletListAdressATLJReacFieldGenTab::setPotentialAT2)
        .def("setPotentialCG", &VerletListAdressATLJReacFieldGenTab::setPotentialCG);
    ;

    class_<VerletListAdressATLJReacFieldGenHarmonic, bases<Interaction> >(
        "interaction_VerletListAdressATLJReacFieldGenHarmonic",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListAdressATLJReacFieldGenHarmonic::setPotentialAT1)
        .def("setPotentialAT2", &VerletListAdressATLJReacFieldGenHarmonic::setPotentialAT2)
        .def("setPotentialCG", &VerletListAdressATLJReacFieldGenHarmonic::setPotentialCG);
    ;

    class_<VerletListAdressCGLennardJones, bases<Interaction> >(
        "interaction_VerletListAdressCGLennardJones",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListAdressCGLennardJones::getVerletList)
        .def("setPotential", &VerletListAdressCGLennardJones::setPotential)
        .def("getPotential", &VerletListAdressCGLennardJones::getPotentialPtr);

    class_<VerletListAdressLennardJones, bases<Interaction> >(
        "interaction_VerletListAdressLennardJones",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJones::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJones::setPotentialCG);
    ;

    class_<VerletListAdressLennardJones2, bases<Interaction> >(
        "interaction_VerletListAdressLennardJones2",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJones2::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJones2::setPotentialCG);
    ;

    class_<VerletListAdressLennardJonesHarmonic, bases<Interaction> >(
        "interaction_VerletListAdressLennardJonesHarmonic",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesHarmonic::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesHarmonic::setPotentialCG);
    ;

    class_<VerletListHadressATLennardJones, bases<Interaction> >(
        "interaction_VerletListHadressATLennardJones",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressATLennardJones::getVerletList)
        .def("setPotential", &VerletListHadressATLennardJones::setPotential)
        .def("getPotential", &VerletListHadressATLennardJones::getPotentialPtr);

    class_<VerletListHadressATLenJonesReacFieldGen, bases<Interaction> >(
        "interaction_VerletListHadressATLenJonesReacFieldGen",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotential1", &VerletListHadressATLenJonesReacFieldGen::setPotential1)
        .def("setPotential2", &VerletListHadressATLenJonesReacFieldGen::setPotential2);
    ;

    class_<VerletListHadressATLJReacFieldGenTab, bases<Interaction> >(
        "interaction_VerletListHadressATLJReacFieldGenTab",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListHadressATLJReacFieldGenTab::setPotentialAT1)
        .def("setPotentialAT2", &VerletListHadressATLJReacFieldGenTab::setPotentialAT2)
        .def("setPotentialCG", &VerletListHadressATLJReacFieldGenTab::setPotentialCG);
    ;

    class_<VerletListHadressATLJReacFieldGenHarmonic, bases<Interaction> >(
        "interaction_VerletListHadressATLJReacFieldGenHarmonic",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT1", &VerletListHadressATLJReacFieldGenHarmonic::setPotentialAT1)
        .def("setPotentialAT2", &VerletListHadressATLJReacFieldGenHarmonic::setPotentialAT2)
        .def("setPotentialCG", &VerletListHadressATLJReacFieldGenHarmonic::setPotentialCG);
    ;

    class_<VerletListHadressCGLennardJones, bases<Interaction> >(
        "interaction_VerletListHadressCGLennardJones",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("getVerletList", &VerletListHadressCGLennardJones::getVerletList)
        .def("setPotential", &VerletListHadressCGLennardJones::setPotential)
        .def("getPotential", &VerletListHadressCGLennardJones::getPotentialPtr);

    class_<VerletListHadressLennardJones, bases<Interaction> >(
        "interaction_VerletListHadressLennardJones",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJones::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJones::setPotentialCG);
    ;

    class_<VerletListHadressLennardJones2, bases<Interaction> >(
        "interaction_VerletListHadressLennardJones2",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJones2::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJones2::setPotentialCG);
    ;

    class_<VerletListHadressLennardJonesHarmonic, bases<Interaction> >(
        "interaction_VerletListHadressLennardJonesHarmonic",
        init<std::shared_ptr<VerletListAdress>, std::shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJonesHarmonic::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJonesHarmonic::setPotentialCG);
    ;

    class_<VerletListHybridLennardJones, bases<Interaction> >(
        "interaction_VerletListHybridLennardJones",
        init<std::shared_ptr<VerletList>, bool>())
        .def("getVerletList", &VerletListHybridLennardJones::getVerletList)
        .def("setPotential", &VerletListHybridLennardJones::setPotential)
        .def("getPotential", &VerletListHybridLennardJones::getPotentialPtr)
        .add_property("scale_factor", &VerletListHybridLennardJones::scaleFactor,
                      &VerletListHybridLennardJones::setScaleFactor)
        .add_property("max_force", &VerletListHybridLennardJones::maxForce,
                      &VerletListHybridLennardJones::setMaxForce);

    class_<CellListLennardJones, bases<Interaction> >("interaction_CellListLennardJones",
                                                      init<std::shared_ptr<storage::Storage> >())
        .def("setPotential", &CellListLennardJones::setPotential);
    ;

    class_<FixedPairListLennardJones, bases<Interaction> >(
        "interaction_FixedPairListLennardJones",
        init<std::shared_ptr<System>, std::shared_ptr<FixedPairList>,
             std::shared_ptr<LennardJones> >())
        .def(init<std::shared_ptr<System>, std::shared_ptr<FixedPairListAdress>,
                  std::shared_ptr<LennardJones> >())
        .def("setPotential", &FixedPairListLennardJones::setPotential)
        .def("getPotential", &FixedPairListLennardJones::getPotential)
        .def("setFixedPairList", &FixedPairListLennardJones::setFixedPairList)
        .def("getFixedPairList", &FixedPairListLennardJones::getFixedPairList);

    class_<FixedPairListTypesLennardJones, bases<Interaction> >(
        "interaction_FixedPairListTypesLennardJones",
        init<std::shared_ptr<System>, std::shared_ptr<FixedPairList> >())
        .def(init<std::shared_ptr<System>, std::shared_ptr<FixedPairListAdress> >())
        .def("setFixedPairList", &FixedPairListTypesLennardJones::setFixedPairList)
        .def("getFixedPairList", &FixedPairListTypesLennardJones::getFixedPairList)
        .def("setPotential", &FixedPairListTypesLennardJones::setPotential)
        .def("getPotential", &FixedPairListTypesLennardJones::getPotentialPtr);

    class_<FixedPairListAdressLennardJones, bases<Interaction> >(
        "interaction_FixedPairListAdressLennardJones",
        init<std::shared_ptr<System>, std::shared_ptr<FixedPairList>, std::shared_ptr<LennardJones>, bool>())
        .def(init<std::shared_ptr<System>, std::shared_ptr<FixedPairListAdress>, std::shared_ptr<LennardJones>,
             bool>())
        .def("setPotential", &FixedPairListAdressLennardJones::setPotential)
        .def("getPotential", &FixedPairListAdressLennardJones::getPotential)
        .def("setFixedPairList", &FixedPairListAdressLennardJones::setFixedPairList)
        .def("getFixedPairList", &FixedPairListAdressLennardJones::getFixedPairList)
        ;
}

}  // namespace interaction
}  // namespace espressopp
