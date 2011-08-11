#include "python.hpp"
#include "Potential.hpp"
#include "logging.hpp"

namespace espresso {
  namespace interaction {

    LOG4ESPP_LOGGER(Potential::theLogger, "Potential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Potential::registerPython() {
        using namespace espresso::python;
        
        real (Potential::*computeEnergy1)(const Real3D& dist) const =
            &Potential::computeEnergy;
            
        real (Potential::*computeEnergy2)(real dist) const =
            &Potential::computeEnergy;
            
        Real3D (Potential::*computeForce)(const Real3D& dist) const =
            &Potential::computeForce;
        
        class_< Potential, boost::noncopyable >
            ("interaction_Potential", no_init)
            .add_property("cutoff",  
                &Potential::getCutoff,
                &Potential::setCutoff)
            .add_property("shift",  
                &Potential::getShift,  
                &Potential::setShift)
            .def("setAutoShift", pure_virtual(&Potential::setAutoShift))
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce))
        ;
    }
  }
}


//     class PythonPotential 
//       : public espresso::python::wrapper< Potential >,
// 	public PotentialBase< PythonPotential > {
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
//     Potential::registerPython() {
//       using namespace espresso::python;
      
//       class_< espresso::potential::PythonPotential, boost::noncopyable >
// 	("interaction_PythonPotential")
// 	.def("getCutoffSqr", pure_virtual(&Potential::getCutoffSqr))
// 	.def("computeEnergy", pure_virtual(&Potential::computeEnergy))
// 	.def("computeForce", pure_virtual(&Potential::computeForce))
// 	;
//     }
