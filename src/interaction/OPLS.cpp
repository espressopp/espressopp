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
#include "OPLS.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    OPLS::registerPython() {
      using namespace espressopp::python;

      class_ <OPLS, bases <DihedralPotential> >
    	("interaction_OPLS", init< real, real, real, real >())
	.add_property("K1", &OPLS::getK1, &OPLS::setK1)
	.add_property("K2", &OPLS::getK2, &OPLS::setK2)
	.add_property("K3", &OPLS::getK3, &OPLS::setK3)
	.add_property("K4", &OPLS::getK4, &OPLS::setK4)
        //set all K at once
    	;

      typedef class FixedQuadrupleListInteractionTemplate <OPLS>
        FixedQuadrupleListOPLS;
      class_ <FixedQuadrupleListOPLS, bases <Interaction> >
        ("interaction_FixedQuadrupleListOPLS",
                init< shared_ptr<System>,
                      shared_ptr<FixedQuadrupleList>,
                      shared_ptr<OPLS> >())
        .def("setPotential", &FixedQuadrupleListOPLS::setPotential);
      ;
    }
  }
}
