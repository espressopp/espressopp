#include "python.hpp"
#include "DihedralPotential.hpp"

namespace espresso {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralPotential::registerPython() {
        using namespace espresso::python;
        
        real (DihedralPotential::*computeEnergy1)
            (const Real3D& dist21, const Real3D& dist32, const Real3D& dist43) const =
                &DihedralPotential::computeEnergy;
        
        real (DihedralPotential::*computeEnergy2)
            (real phi) const =
                &DihedralPotential::computeEnergy;
        
        void (DihedralPotential::*computeForce1)
            (Real3D& force1, Real3D& force2, Real3D& force3, Real3D& force4,
            const Real3D& dist21,const Real3D& dist32, const Real3D& dist43) const =
                &DihedralPotential::computeForce;
        
        real (DihedralPotential::*computeForce2)
            (real phi) const =
                &DihedralPotential::computeForce;

        class_< DihedralPotential, boost::noncopyable >
            ("interaction_DihedralPotential", no_init)
            .add_property("cutoff",
                &DihedralPotential::getCutoff,
                &DihedralPotential::setCutoff)
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce1))
            .def("computeForce", pure_virtual(computeForce2))
        ;
    }
  }
}
