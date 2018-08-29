/*
  Copyright (C) 2012,2013,2018
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

    real
    DihedralPotential::computePhi(const Real3D& r21,
                 const Real3D& r32,
                 const Real3D& r43) {
      Real3D rijjk = r21.cross(r32); // [r21 x r32]
      Real3D rjkkn = r32.cross(r43); // [r32 x r43]

      real rijjk_abs = rijjk.abs();
      real rjkkn_abs = rjkkn.abs();

      real inv_rijjk = 1.0 / rijjk_abs;
      real inv_rjkkn = 1.0 / rjkkn_abs;

      // cosine between planes
      real cos_phi = (rijjk * rjkkn) * (inv_rijjk * inv_rjkkn);
      if (cos_phi > 1.0) cos_phi = 1.0;
      else if (cos_phi < -1.0) cos_phi = -1.0;

      real phi = acos(cos_phi);
      //get sign of phi
      //positive if (rij x rjk) x (rjk x rkn) is in the same direction as rjk, negative otherwise (see DLPOLY manual)
      Real3D rcross = rijjk.cross(rjkkn); //(rij x rjk) x (rjk x rkn)
      real signcheck = rcross * r32;
      if (signcheck < 0.0) phi *= -1.0;

      return phi;
    }

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
            .add_property("colVarBondList",
                &DihedralPotential::getColVarBondList,
                &DihedralPotential::setColVarBondList)
            .add_property("colVarAngleList",
                &DihedralPotential::getColVarAngleList,
                &DihedralPotential::setColVarAngleList)
            .add_property("colVar",
                &DihedralPotential::getColVar)
            .def("computeEnergy", pure_virtual(computeEnergy1))
            .def("computeEnergy", pure_virtual(computeEnergy2))
            .def("computeForce", pure_virtual(computeForce1))
            .def("computeForce", pure_virtual(computeForce2))
            .def("computePhi", &DihedralPotential::computePhi)
        ;
    }
  }
}
