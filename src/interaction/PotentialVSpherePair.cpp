#include "python.hpp"
#include "PotentialVSpherePair.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(PotentialVSpherePair::theLogger, "PotentialVSpherePair");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void PotentialVSpherePair::registerPython() {
        using namespace espresso::python;
        
        real (PotentialVSpherePair::*computeEnergy1)(const Real3D& dist, real& sigmaij) const =
            &PotentialVSpherePair::computeEnergy;
            
        real (PotentialVSpherePair::*computeEnergy2)(real dist, real sigmaij) const =
            &PotentialVSpherePair::computeEnergy;
            
        python::list (PotentialVSpherePair::*computeForce)(const Real3D& dist, const real& sigmai, const real& sigmaj) const =
            &PotentialVSpherePair::computeForce;
        
        class_< PotentialVSpherePair, boost::noncopyable >
            ("interaction_PotentialVSpherePair", no_init)
            .add_property("cutoff",  
                &PotentialVSpherePair::getCutoff,
                &PotentialVSpherePair::setCutoff)
            .add_property("shift",  
                &PotentialVSpherePair::getShift,
                &PotentialVSpherePair::setShift)
            .def("setAutoShift", pure_virtual(&PotentialVSpherePair::setAutoShift))
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce))
        ;
    }
  }
}
