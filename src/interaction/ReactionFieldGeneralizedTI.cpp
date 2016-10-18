/*
  Copyright (C) 2012,2013,2016
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
#include "ReactionFieldGeneralizedTI.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {

        //typedef class VerletListInteractionTemplate<ReactionFieldGeneralizedTI>
         //   VerletListReactionFieldGeneralizedTI;

        typedef class VerletListAdressInteractionTemplate<ReactionFieldGeneralizedTI, Tabulated>
            VerletListAdressReactionFieldGeneralizedTI;

        //typedef class VerletListHadressInteractionTemplate<ReactionFieldGeneralizedTI, Tabulated>
        //    VerletListHadressReactionFieldGeneralizedTI;
        
        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void ReactionFieldGeneralizedTI::registerPython() {
            using namespace espressopp::python;

            void (ReactionFieldGeneralizedTI::*pyAddPid)(longint pid) = &ReactionFieldGeneralizedTI::addPid;

            class_<ReactionFieldGeneralizedTI, bases<Potential> >
                ("interaction_ReactionFieldGeneralizedTI", init< real, real, real, real, real, real, bool>())
                .def_pickle(ReactionFieldGeneralizedTI_pickle())
                .def("addPid",pyAddPid)
                .add_property("prefactor", &ReactionFieldGeneralizedTI::getPrefactor, &ReactionFieldGeneralizedTI::setPrefactor)
            ;

            //class_<VerletListReactionFieldGeneralizedTI, bases<Interaction> >
            //    ("interaction_VerletListReactionFieldGeneralizedTI", init< shared_ptr<VerletList> >())
            //    .def("setPotential", &VerletListReactionFieldGeneralizedTI::setPotential)
            //    .def("getPotential", &VerletListReactionFieldGeneralizedTI::getPotentialPtr)
            //;

            class_<VerletListAdressReactionFieldGeneralizedTI, bases<Interaction> >
                ("interaction_VerletListAdressReactionFieldGeneralizedTI",
                        init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
                .def("setPotentialAT", &VerletListAdressReactionFieldGeneralizedTI::setPotentialAT)
                .def("setPotentialCG", &VerletListAdressReactionFieldGeneralizedTI::setPotentialCG);
            ;

            //class_<VerletListHadressReactionFieldGeneralizedTI, bases<Interaction> >
            //    ("interaction_VerletListHadressReactionFieldGeneralizedTI",
            //            init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
            //    .def("setPotentialAT", &VerletListHadressReactionFieldGeneralizedTI::setPotentialAT)
            //    .def("setPotentialCG", &VerletListHadressReactionFieldGeneralizedTI::setPotentialCG);
            //;
            
        }
    }
}
