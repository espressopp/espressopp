/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"
#include "Potential.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(Potential::theLogger, "Potential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Potential::registerPython() {
        using namespace espressopp::python;
        
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
//       : public espressopp::python::wrapper< Potential >,
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
//       using namespace espressopp::python;
      
//       class_< espressopp::potential::PythonPotential, boost::noncopyable >
// 	("interaction_PythonPotential")
// 	.def("getCutoffSqr", pure_virtual(&Potential::getCutoffSqr))
// 	.def("computeEnergy", pure_virtual(&Potential::computeEnergy))
// 	.def("computeForce", pure_virtual(&Potential::computeForce))
// 	;
//     }
