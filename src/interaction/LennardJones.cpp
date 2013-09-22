#include "python.hpp"
#include "LennardJones.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJones>
        VerletListLennardJones;
    typedef class VerletListAdressInteractionTemplate <LennardJones, Tabulated>
        VerletListAdressLennardJones;
    typedef class VerletListAdressInteractionTemplate <LennardJones, LennardJones>
        VerletListAdressLennardJones2;
    typedef class VerletListHadressInteractionTemplate <LennardJones, Tabulated>
        VerletListHadressLennardJones;
    typedef class VerletListHadressInteractionTemplate <LennardJones, LennardJones>
        VerletListHadressLennardJones2;
    typedef class CellListAllPairsInteractionTemplate <LennardJones> 
        CellListLennardJones;
    typedef class FixedPairListInteractionTemplate <LennardJones> 
        FixedPairListLennardJones;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJones::registerPython() {
      using namespace espresso::python;

      class_< LennardJones, bases< Potential > >
    	("interaction_LennardJones", init< real, real, real >())
	    .def(init< real, real, real, real >())
    	.add_property("sigma", &LennardJones::getSigma, &LennardJones::setSigma)
    	.add_property("epsilon", &LennardJones::getEpsilon, &LennardJones::setEpsilon)
        .def_pickle(LennardJones_pickle())

      ;

      class_< VerletListLennardJones, bases< Interaction > > 
        ("interaction_VerletListLennardJones", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJones::getVerletList)
        .def("setPotential", &VerletListLennardJones::setPotential)
        .def("getPotential", &VerletListLennardJones::getPotentialPtr)
      ;

      class_< VerletListAdressLennardJones, bases< Interaction > >
        ("interaction_VerletListAdressLennardJones",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJones::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJones::setPotentialCG);
      ;
      
      class_< VerletListAdressLennardJones2, bases< Interaction > >
        ("interaction_VerletListAdressLennardJones2",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJones2::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJones2::setPotentialCG);
      ;

      class_< VerletListHadressLennardJones, bases< Interaction > >
        ("interaction_VerletListHadressLennardJones",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress>, bool >())
        .def("setPotentialAT", &VerletListHadressLennardJones::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJones::setPotentialCG);
      ;
      
      class_< VerletListHadressLennardJones2, bases< Interaction > >
        ("interaction_VerletListHadressLennardJones2",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress>, bool >())
        .def("setPotentialAT", &VerletListHadressLennardJones2::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJones2::setPotentialCG);
      ;
      
      class_< CellListLennardJones, bases< Interaction > > 
        ("interaction_CellListLennardJones", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJones::setPotential);
	  ;

      class_< FixedPairListLennardJones, bases< Interaction > >
        ("interaction_FixedPairListLennardJones",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJones> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJones> >())
          .def("setPotential", &FixedPairListLennardJones::setPotential)
          .def("setFixedPairList", &FixedPairListLennardJones::setFixedPairList)
          .def("getFixedPairList", &FixedPairListLennardJones::getFixedPairList)
      ;
    }
    
  }
}
