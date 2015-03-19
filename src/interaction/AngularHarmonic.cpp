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
#include "AngularHarmonic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularHarmonic::registerPython() {
      using namespace espressopp::python;

      class_< AngularHarmonic, bases< AngularPotential > >
    	("interaction_AngularHarmonic", init< real, real >())
	.add_property("K", &AngularHarmonic::getK, &AngularHarmonic::setK)
	.add_property("theta0", &AngularHarmonic::getTheta0, &AngularHarmonic::setTheta0)
    	;

      typedef class FixedTripleListInteractionTemplate<AngularHarmonic>
        FixedTripleListAngularHarmonic;
        
      class_ <FixedTripleListAngularHarmonic, bases <Interaction> >
        ("interaction_FixedTripleListAngularHarmonic",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleList>,
                shared_ptr<AngularHarmonic> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedTripleListAdress>, shared_ptr<AngularHarmonic> >())
        .def("setPotential", &FixedTripleListAngularHarmonic::setPotential)
        .def("getFixedTripleList", &FixedTripleListAngularHarmonic::getFixedTripleList);
      ;
    }
  }
}
