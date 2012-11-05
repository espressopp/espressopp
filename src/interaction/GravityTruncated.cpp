#include "python.hpp"
#include "GravityTruncated.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <GravityTruncated> VerletListGravityTruncated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void GravityTruncated::registerPython() {
      using namespace espresso::python;

      class_< GravityTruncated, bases< Potential > >
        ("interaction_GravityTruncated", init< >())
        .def(init< real, real >())
        .add_property("prefactor", &GravityTruncated::getPrefactor, &GravityTruncated::setPrefactor)
      ;

      class_< VerletListGravityTruncated, bases< Interaction > >
        ("interaction_VerletListGravityTruncated", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListGravityTruncated::getVerletList)
        .def("setPotential", &VerletListGravityTruncated::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListGravityTruncated::getPotential, return_value_policy< reference_existing_object >())
      ;
    }
  }
}
