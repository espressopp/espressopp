#include "python.hpp"
#include "StillingerWeberTripleTerm.hpp"

//#include "VerletListTripleInteractionTemplate.hpp"
//#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    typedef class VerletListTripleInteractionTemplate <StillingerWeberTripleTerm>
        VerletListStillingerWeberTripleTerm;
    typedef class FixedTripleListInteractionTemplate <StillingerWeberTripleTerm>
        FixedTripleListStillingerWeberTripleTerm;
    
    void 
    StillingerWeberTripleTerm::registerPython() {
      using namespace espresso::python;

      class_< StillingerWeberTripleTerm, bases< AngularPotential > >(
            "interaction_StillingerWeberTripleTerm",
            init< real, real, real, real, real, real, real, real, real >()
      )
      .add_property("A", &StillingerWeberTripleTerm::getA, &StillingerWeberTripleTerm::setA)
      .add_property("B", &StillingerWeberTripleTerm::getB, &StillingerWeberTripleTerm::setB)
      .add_property("p", &StillingerWeberTripleTerm::getP, &StillingerWeberTripleTerm::setP)
      .add_property("q", &StillingerWeberTripleTerm::getQ, &StillingerWeberTripleTerm::setQ)
      .add_property("gamma", &StillingerWeberTripleTerm::getGamma, &StillingerWeberTripleTerm::setGamma)
      .add_property("theta0", &StillingerWeberTripleTerm::getTheta0, &StillingerWeberTripleTerm::setTheta0)
      .add_property("lambda", &StillingerWeberTripleTerm::getLambda, &StillingerWeberTripleTerm::setLambda)
      .add_property("epsilon", &StillingerWeberTripleTerm::getEpsilon, &StillingerWeberTripleTerm::setEpsilon)
      .add_property("sigma", &StillingerWeberTripleTerm::getSigma, &StillingerWeberTripleTerm::setSigma)
      ;

      
      class_< VerletListStillingerWeberTripleTerm, bases< Interaction > >(
        "interaction_VerletListStillingerWeberTripleTerm",
        init< shared_ptr<System>, shared_ptr<VerletListTriple>, shared_ptr <StillingerWeberTripleTerm> >()
      )
      .def("getVerletListTriple", &VerletListStillingerWeberTripleTerm::getVerletListTriple)
      .def("setPotential", &VerletListStillingerWeberTripleTerm::setPotential)
      //.def("getPotential", &VerletListStillingerWeberTripleTerm::getPotential,
      //        return_value_policy< reference_existing_object >())
      ;

      class_< FixedTripleListStillingerWeberTripleTerm, bases< Interaction > >(
          "interaction_FixedTripleListStillingerWeberTripleTerm",
          init<shared_ptr<System>, shared_ptr<FixedTripleList>, shared_ptr<StillingerWeberTripleTerm> >()
      )
      .def("setPotential", &FixedTripleListStillingerWeberTripleTerm::setPotential)
      .def("getFixedTripleList", &FixedTripleListStillingerWeberTripleTerm::getFixedTripleList)
      ;
    }
  }
}
