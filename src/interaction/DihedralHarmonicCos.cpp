#include "python.hpp"
#include "DihedralHarmonicCos.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    DihedralHarmonicCos::registerPython() {
      using namespace espresso::python;

      class_ <DihedralHarmonicCos, bases <DihedralPotential> >
    	("interaction_DihedralHarmonicCos", init< real, real >())
        .add_property("K", &DihedralHarmonicCos::getK, &DihedralHarmonicCos::setK)
        .add_property("phi0", &DihedralHarmonicCos::getPhi0, &DihedralHarmonicCos::setPhi0)
    	;

      typedef class FixedQuadrupleListInteractionTemplate <DihedralHarmonicCos>
      FixedQuadrupleListDihedralHarmonicCos;
      class_ <FixedQuadrupleListDihedralHarmonicCos, bases <Interaction> >
        ("interaction_FixedQuadrupleListDihedralHarmonicCos",
                  init< shared_ptr<System>,
                        shared_ptr<FixedQuadrupleList>,
                        shared_ptr<DihedralHarmonicCos> >())
        .def("setPotential", &FixedQuadrupleListDihedralHarmonicCos::setPotential)
        .def("getFixedQuadrupleList", &FixedQuadrupleListDihedralHarmonicCos::getFixedQuadrupleList)
        ;
    }
  }
}
