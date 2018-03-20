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

#include "DihedralHarmonicUniqueCos.hpp"
#include "FixedQuadrupleAngleListInteractionTemplate.hpp"
#include "python.hpp"

namespace espressopp {
namespace interaction {
//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void DihedralHarmonicUniqueCos::registerPython() {
  using namespace espressopp::python;

  class_<DihedralHarmonicUniqueCos, bases<DihedralUniquePotential> >(
      "interaction_DihedralHarmonicUniqueCos", init<real>())
      .add_property("K", &DihedralHarmonicUniqueCos::getK, &DihedralHarmonicUniqueCos::setK);

  typedef class FixedQuadrupleAngleListInteractionTemplate<DihedralHarmonicUniqueCos>
      FixedQuadrupleAngleListDihedralHarmonicUniqueCos;
  class_<FixedQuadrupleAngleListDihedralHarmonicUniqueCos, bases<Interaction> >(
      "interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos",
      init<shared_ptr<System>, shared_ptr<FixedQuadrupleAngleList>,
           shared_ptr<DihedralHarmonicUniqueCos> >())
      .def("setPotential", &FixedQuadrupleAngleListDihedralHarmonicUniqueCos::setPotential)
      .def("getFixedQuadrupleAngleList",
           &FixedQuadrupleAngleListDihedralHarmonicUniqueCos::getFixedQuadrupleAngleList);
}
}
}
