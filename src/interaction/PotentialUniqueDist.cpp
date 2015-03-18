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
#include "PotentialUniqueDist.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(PotentialUniqueDist::theLogger, "PotentialUniqueDist");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void PotentialUniqueDist::registerPython() {
        using namespace espressopp::python;
        
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
//       : public espressopp::python::wrapper< PotentialUniqueDist >,
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
//       using namespace espressopp::python;
      
//       class_< espressopp::PotentialUniqueDist::PythonPotentialUniqueDist, boost::noncopyable >
// 	("interaction_PythonPotentialUniqueDist")
// 	.def("getCutoffSqr", pure_virtual(&PotentialUniqueDist::getCutoffSqr))
// 	.def("computeEnergy", pure_virtual(&PotentialUniqueDist::computeEnergy))
// 	.def("computeForce", pure_virtual(&PotentialUniqueDist::computeForce))
// 	;
//     }
