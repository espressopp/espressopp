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
#include "DihedralPotential.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(DihedralPotential::theLogger, "DihedralPotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralPotential::registerPython() {
        using namespace espressopp::python;
        
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
