#include "python.hpp"
#include "DihedralUniquePotential.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(DihedralUniquePotential::theLogger, "DihedralUniquePotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralUniquePotential::registerPython() {
        using namespace espresso::python;
        
        real (DihedralUniquePotential::*computeEnergy1)
            (const Real3D& r21, const Real3D& r32, const Real3D& r43, const real phi0) const =
                &DihedralUniquePotential::computeEnergy;
        
        real (DihedralUniquePotential::*computeEnergy2)
            (real phi, real phi0) const =
                &DihedralUniquePotential::computeEnergy;
        
        void (DihedralUniquePotential::*computeForce1)
            (Real3D& force1, Real3D& force2, Real3D& force3, Real3D& force4,
            const Real3D& r21,const Real3D& r32, const Real3D& r43, const real phi0) const =
                &DihedralUniquePotential::computeForce;
        
        real (DihedralUniquePotential::*computeForce2)
            (real phi, real phi0) const =
                &DihedralUniquePotential::computeForce;

        class_< DihedralUniquePotential, boost::noncopyable >
            ("interaction_DihedralUniquePotential", no_init)
            .add_property("cutoff",
                &DihedralUniquePotential::getCutoff,
                &DihedralUniquePotential::setCutoff)
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce1))
            .def("computeForce", pure_virtual(computeForce2))
        ;
    }
  }
}
