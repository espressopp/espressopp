#include "python.hpp"
#include "AngularPotential.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(AngularPotential::theLogger, "AngularPotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void AngularPotential::registerPython() {
        using namespace espresso::python;
        
        real (AngularPotential::*computeEnergy1)
            (const Real3D& dist12, const Real3D& dist32) const =
                &AngularPotential::computeEnergy;
        
        real (AngularPotential::*computeEnergy2)
            (real theta) const =
                &AngularPotential::computeEnergy;
        
        void (AngularPotential::*computeForce1)
            (Real3D& force12, Real3D& force32,
            const Real3D& dist12, const Real3D& dist32) const =
                &AngularPotential::computeForce;
        
        real (AngularPotential::*computeForce2)(real theta) const =
                &AngularPotential::computeForce;
        
        class_< AngularPotential, boost::noncopyable >
            ("interaction_AngularPotential", no_init)
            .add_property("cutoff",
                &AngularPotential::getCutoff,
                &AngularPotential::setCutoff)
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce1))
            .def("computeForce", pure_virtual(computeForce2))
        ;
    }
  }
}
