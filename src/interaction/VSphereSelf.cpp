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
#include "VSphereSelf.hpp"
#include "VSphereSelfInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
      typedef class VSphereSelfInteractionTemplate< VSphereSelf >
      SelfVSphere;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    VSphereSelf::registerPython() {
      using namespace espressopp::python;

      class_< VSphereSelf, bases< Potential > >("interaction_VSphereSelf",
      init< real, real, real, int, real >())
	  .def(init< real, real, real, int, real, real >())
	  .add_property("e1", &VSphereSelf::gete1, &VSphereSelf::sete1)
	  .add_property("a1", &VSphereSelf::geta1, &VSphereSelf::seta1)
	  .add_property("a2", &VSphereSelf::geta2, &VSphereSelf::seta2)
	  .add_property("Nb", &VSphereSelf::getNb, &VSphereSelf::setNb)
      ;

      class_< SelfVSphere, bases< Interaction > >("interaction_SelfVSphere",
      init< shared_ptr<System>, shared_ptr<VSphereSelf> >())
      .def("setPotential", &SelfVSphere::setPotential)
      .def("getPotential", &SelfVSphere::getPotential)
      ;
    }

  }
}
