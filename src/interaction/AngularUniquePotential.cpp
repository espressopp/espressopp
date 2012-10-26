#include "python.hpp"
#include "AngularUniquePotential.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(AngularUniquePotential::theLogger, "AngularUniquePotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void AngularUniquePotential::registerPython() {
      using namespace espresso::python;

      real (AngularUniquePotential::*computeEnergy1)
          (const Real3D& dist12, const Real3D& dist32, real cos0) const =
              &AngularUniquePotential::computeEnergy;

      real (AngularUniquePotential::*computeEnergy2)
          (real theta, real cos0) const =
              &AngularUniquePotential::computeEnergy;

      void (AngularUniquePotential::*computeForce1)
          (Real3D& force12, Real3D& force32,
          const Real3D& dist12, const Real3D& dist32, real cos0) const =
              &AngularUniquePotential::computeForce;

      real (AngularUniquePotential::*computeForce2)(real theta, real cos0) const =
              &AngularUniquePotential::computeForce;

      class_< AngularUniquePotential, boost::noncopyable >
          ("interaction_AngularUniquePotential", no_init)
          .add_property("cutoff",
              &AngularUniquePotential::getCutoff,
              &AngularUniquePotential::setCutoff)
          .def("computeEnergy", pure_virtual(computeEnergy1))
          .def("computeEnergy", pure_virtual(computeEnergy2))
          .def("computeForce", pure_virtual(computeForce1))
          .def("computeForce", pure_virtual(computeForce2))
      ;
    }
  }
}
