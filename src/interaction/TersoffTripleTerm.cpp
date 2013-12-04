#include "python.hpp"
#include "TersoffTripleTerm.hpp"

#include "VerletListTripleInteractionTemplate.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    typedef class VerletListTripleInteractionTemplate <TersoffTripleTerm>
        VerletListTersoffTripleTerm;
    typedef class FixedTripleListInteractionTemplate <TersoffTripleTerm>
        FixedListTersoffTripleTerm;
    
    void 
    TersoffTripleTerm::registerPython() {
      using namespace espresso::python;

      class_< TersoffTripleTerm, bases< AngularPotential > >(
            "interaction_TersoffTripleTerm",
            init< real, real, real, real, real, real, real, real, real, real, real, real, real, real >()
      )
      .add_property("B", &TersoffTripleTerm::getB,
              &TersoffTripleTerm::setB)
      .add_property("lambda2", &TersoffTripleTerm::getLambda2,
              &TersoffTripleTerm::setLambda2)
      .add_property("R", &TersoffTripleTerm::getR,
              &TersoffTripleTerm::setR)
      .add_property("D", &TersoffTripleTerm::getD,
              &TersoffTripleTerm::setD)
      .add_property("n", &TersoffTripleTerm::getN,
              &TersoffTripleTerm::setN)
      .add_property("beta", &TersoffTripleTerm::getBeta,
              &TersoffTripleTerm::setBeta)
      .add_property("m", &TersoffTripleTerm::getM,
              &TersoffTripleTerm::setM)
      .add_property("lambda3", &TersoffTripleTerm::getLambda3,
              &TersoffTripleTerm::setLambda3)
      .add_property("gamma", &TersoffTripleTerm::getGamma,
              &TersoffTripleTerm::setGamma)
      .add_property("c", &TersoffTripleTerm::getC,
              &TersoffTripleTerm::setC)
      .add_property("d", &TersoffTripleTerm::getd,
              &TersoffTripleTerm::setd)
      .add_property("theta0", &TersoffTripleTerm::getTheta0,
              &TersoffTripleTerm::setTheta0)
      .add_property("cutoff1", &TersoffTripleTerm::getCutoff1,
              &TersoffTripleTerm::setCutoff1)
      .add_property("cutoff2", &TersoffTripleTerm::getCutoff2,
              &TersoffTripleTerm::setCutoff2)
      ;

      
      class_< VerletListTersoffTripleTerm, bases< Interaction > >(
        "interaction_VerletListTersoffTripleTerm",
        init< shared_ptr<System>, 
              shared_ptr<VerletListTriple> >()
      )
      .def("getVerletListTriple", &VerletListTersoffTripleTerm::getVerletListTriple)
      .def("setPotential", &VerletListTersoffTripleTerm::setPotential,
              return_value_policy< reference_existing_object >())
      .def("getPotential", &VerletListTersoffTripleTerm::getPotential,
              return_value_policy< reference_existing_object >())
      ;

      class_< FixedListTersoffTripleTerm, bases< Interaction > >(
          "interaction_FixedTripleListTersoffTripleTerm",
          init<shared_ptr<System>, shared_ptr<FixedTripleList>, shared_ptr<TersoffTripleTerm> >()
      )
      .def("setPotential", &FixedListTersoffTripleTerm::setPotential)
      .def("getFixedTripleList", &FixedListTersoffTripleTerm::getFixedTripleList)
      ;
    }
  }
}
