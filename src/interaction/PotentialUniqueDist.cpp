#include "python.hpp"
#include "PotentialUniqueDist.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(PotentialUniqueDist::theLogger, "PotentialUniqueDist");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void PotentialUniqueDist::registerPython() {
        using namespace espresso::python;
        
        real (PotentialUniqueDist::*computeEnergy1)(const Real3D& dist, const real curDist) const =
            &PotentialUniqueDist::computeEnergy;
            
        real (PotentialUniqueDist::*computeEnergy2)(real dist, real curDist) const =
            &PotentialUniqueDist::computeEnergy;
            
        Real3D (PotentialUniqueDist::*computeForce)(const Real3D& dist, const real curDist) const =
            &PotentialUniqueDist::computeForce;
        
        class_< PotentialUniqueDist, boost::noncopyable >
            ("interaction_PotentialUniqueDist", no_init)
            .add_property("cutoff",  
                &PotentialUniqueDist::getCutoff,
                &PotentialUniqueDist::setCutoff)
            .add_property("shift",  
                &PotentialUniqueDist::getShift,  
                &PotentialUniqueDist::setShift)
            .def("setAutoShift", pure_virtual(&PotentialUniqueDist::setAutoShift))
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce))
        ;
    }
  }
}


//     class PythonPotentialUniqueDist 
//       : public espresso::python::wrapper< PotentialUniqueDist >,
// 	public PotentialUniqueDistBase< PythonPotentialUniqueDist > {
//     public:
//       real _getCutoffSqr() const {
// 	return get_override("getCutoffSqr")();
//       }
      
//       real _computeEnergy(const Real3D dist) const {
// 	return get_override("computeEnergy")(dist);
//       }
      
//       Real3D _computeForce(const Real3D dist) const {
// 	return get_override("computeForce")(dist);
//       }
//     };

//     void
//     PotentialUniqueDist::registerPython() {
//       using namespace espresso::python;
      
//       class_< espresso::PotentialUniqueDist::PythonPotentialUniqueDist, boost::noncopyable >
// 	("interaction_PythonPotentialUniqueDist")
// 	.def("getCutoffSqr", pure_virtual(&PotentialUniqueDist::getCutoffSqr))
// 	.def("computeEnergy", pure_virtual(&PotentialUniqueDist::computeEnergy))
// 	.def("computeForce", pure_virtual(&PotentialUniqueDist::computeForce))
// 	;
//     }
