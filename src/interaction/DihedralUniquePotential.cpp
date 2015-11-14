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
#include "DihedralUniquePotential.hpp"
#include "logging.hpp"

namespace espressopp {
  namespace interaction {

    LOG4ESPP_LOGGER(DihedralUniquePotential::theLogger, "DihedralUniquePotential");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    DihedralUniquePotential::registerPython() {
        using namespace espressopp::python;
        
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
