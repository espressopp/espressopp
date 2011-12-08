#include "python.hpp"
#include "ReactionFieldGeneralized.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
//#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
    namespace interaction {

        typedef class VerletListInteractionTemplate<ReactionFieldGeneralized>
            VerletListReactionFieldGeneralized;

        typedef class VerletListAdressInteractionTemplate<ReactionFieldGeneralized, Tabulated>
                    VerletListAdressReactionFieldGeneralized;

        typedef class CellListAllPairsInteractionTemplate<ReactionFieldGeneralized>
            CellListReactionFieldGeneralized;

        /*typedef class FixedPairListInteractionTemplate<ReactionFieldGeneralized>
            FixedPairListReactionFieldGeneralized;*/

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void ReactionFieldGeneralized::registerPython() {
            using namespace espresso::python;

            class_<ReactionFieldGeneralized, bases<Potential> >
                ("interaction_ReactionFieldGeneralized", init<real, real, real, real, real>())
                .def(init<real, real, real, real, real, bool>())
                .add_property("qq", &ReactionFieldGeneralized::getQQ, &ReactionFieldGeneralized::setQQ)
            ;

            class_<VerletListReactionFieldGeneralized, bases<Interaction> >
                ("interaction_VerletListReactionFieldGeneralized", init< shared_ptr<VerletList> >())
                .def("setPotential", &VerletListReactionFieldGeneralized::setPotential);
            ;

            class_<VerletListAdressReactionFieldGeneralized, bases<Interaction> >
                ("interaction_VerletListAdressReactionFieldGeneralized",
                        init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleList> >())
                .def("setPotentialAT", &VerletListAdressReactionFieldGeneralized::setPotentialAT)
                .def("setPotentialCG", &VerletListAdressReactionFieldGeneralized::setPotentialCG);
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
