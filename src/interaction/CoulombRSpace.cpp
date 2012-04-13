#include "python.hpp"
#include "CoulombRSpace.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"

// currently just Verlet list

namespace espresso {
  namespace interaction {

    typedef class VerletListInteractionTemplate <CoulombRSpace> VerletListCoulombRSpace;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void CoulombRSpace::registerPython() {
      using namespace espresso::python;

      class_< CoulombRSpace, bases< Potential > >
        ("interaction_CoulombRSpace", init< >())
        .def(init< real, real, real >())
        .add_property("alpha", &CoulombRSpace::getAlpha, &CoulombRSpace::setAlpha)
        .add_property("prefactor", &CoulombRSpace::getPrefactor, &CoulombRSpace::setPrefactor)
      ;

      class_< VerletListCoulombRSpace, bases< Interaction > >
        ("interaction_VerletListCoulombRSpace", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListCoulombRSpace::getVerletList)
        .def("setPotential", &VerletListCoulombRSpace::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListCoulombRSpace::getPotential, return_value_policy< reference_existing_object >())
      ;
    }
    
  }
}
