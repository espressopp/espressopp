#include "python.hpp"
#include "DihedralHarmonicUniqueCos.hpp"
#include "FixedQuadrupleAngleListInteractionTemplate.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    DihedralHarmonicUniqueCos::registerPython() {
      using namespace espresso::python;

      class_ <DihedralHarmonicUniqueCos, bases <DihedralUniquePotential> >
    	("interaction_DihedralHarmonicUniqueCos", init< real >())
        .add_property("K", &DihedralHarmonicUniqueCos::getK, &DihedralHarmonicUniqueCos::setK)
    	;

      typedef class FixedQuadrupleAngleListInteractionTemplate <DihedralHarmonicUniqueCos>
      FixedQuadrupleAngleListDihedralHarmonicUniqueCos;
      class_ <FixedQuadrupleAngleListDihedralHarmonicUniqueCos, bases <Interaction> >
        ("interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos",
                  init< shared_ptr<System>,
                        shared_ptr<FixedQuadrupleAngleList>,
                        shared_ptr<DihedralHarmonicUniqueCos> >())
        .def("setPotential", &FixedQuadrupleAngleListDihedralHarmonicUniqueCos::setPotential)
        .def("getFixedQuadrupleAngleList", 
              &FixedQuadrupleAngleListDihedralHarmonicUniqueCos::getFixedQuadrupleAngleList)
        ;
    }
  }
}
