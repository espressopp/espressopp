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
#include "ReactionFieldGeneralized.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
//#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {

        typedef class VerletListInteractionTemplate<ReactionFieldGeneralized>
            VerletListReactionFieldGeneralized;

        typedef class VerletListAdressInteractionTemplate<ReactionFieldGeneralized, Tabulated>
            VerletListAdressReactionFieldGeneralized;

        typedef class VerletListHadressInteractionTemplate<ReactionFieldGeneralized, Tabulated>
            VerletListHadressReactionFieldGeneralized;
        
        typedef class CellListAllPairsInteractionTemplate<ReactionFieldGeneralized>
            CellListReactionFieldGeneralized;

        /*typedef class FixedPairListInteractionTemplate<ReactionFieldGeneralized>
            FixedPairListReactionFieldGeneralized;*/

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void ReactionFieldGeneralized::registerPython() {
            using namespace espressopp::python;

            class_<ReactionFieldGeneralized, bases<Potential> >
                ("interaction_ReactionFieldGeneralized", init< real, real, real, real, real>())
                .def(init< real, real, real, real, real, real>())
                .def_pickle(ReactionFieldGeneralized_pickle())
                .add_property("prefactor", &ReactionFieldGeneralized::getPrefactor, &ReactionFieldGeneralized::setPrefactor)
            ;

            class_<VerletListReactionFieldGeneralized, bases<Interaction> >
                ("interaction_VerletListReactionFieldGeneralized", init< shared_ptr<VerletList> >())
                .def("setPotential", &VerletListReactionFieldGeneralized::setPotential)
                .def("getPotential", &VerletListReactionFieldGeneralized::getPotentialPtr)
            ;

            class_<VerletListAdressReactionFieldGeneralized, bases<Interaction> >
                ("interaction_VerletListAdressReactionFieldGeneralized",
                        init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
                .def("setPotentialAT", &VerletListAdressReactionFieldGeneralized::setPotentialAT)
                .def("setPotentialCG", &VerletListAdressReactionFieldGeneralized::setPotentialCG);
            ;

            class_<VerletListHadressReactionFieldGeneralized, bases<Interaction> >
                ("interaction_VerletListHadressReactionFieldGeneralized",
                        init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
                .def("setPotentialAT", &VerletListHadressReactionFieldGeneralized::setPotentialAT)
                .def("setPotentialCG", &VerletListHadressReactionFieldGeneralized::setPotentialCG);
            ;
            
            class_<CellListReactionFieldGeneralized, bases<Interaction> >
                ("interaction_CellListReactionFieldGeneralized", init<shared_ptr<storage::Storage> >())
                .def("setPotential", &CellListReactionFieldGeneralized::setPotential);
            ;

            /*
              class_<FixedPairListReactionFieldGeneralized, bases<Interaction> >
                ("interaction_FixedPairListReactionFieldGeneralized",
                  init<shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<ReactionFieldGeneralized> >())
                .def("setPotential", &FixedPairListReactionFieldGeneralized::setPotential);
                ;*/
        }
    }
}
