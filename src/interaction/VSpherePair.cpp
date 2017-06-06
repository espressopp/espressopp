/*
  Copyright (C) 2012,2013,2017
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
#include "VSpherePair.hpp"
#include "Tabulated.hpp"
#include "VerletListVSphereInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListVSphereInteractionTemplate <VSpherePair>
        VerletListVSpherePair;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    VSpherePair::registerPython() {
      using namespace espressopp::python;

      class_< VSpherePair, bases< Potential > >
    	("interaction_VSpherePair", init< real, real >())
	    .def(init< real, real, real >())
    	.add_property("epsilon", &VSpherePair::getEpsilon, &VSpherePair::setEpsilon)
        .def_pickle(VSpherePair_pickle())
      ;

      class_< VerletListVSpherePair, bases< Interaction > >
        ("interaction_VerletListVSpherePair", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListVSpherePair::getVerletList)
        .def("setPotential", &VerletListVSpherePair::setPotential)
        .def("getPotential", &VerletListVSpherePair::getPotentialPtr)
      ;
      
    }
  }
}
