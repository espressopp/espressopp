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
#include "AngularUniqueHarmonic.hpp"
#include "FixedTripleAngleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularUniqueHarmonic::registerPython() {
      using namespace espressopp::python;

      class_< AngularUniqueHarmonic, bases< AngularUniquePotential > >
    	("interaction_AngularUniqueHarmonic", init< real >())
          .add_property("K", &AngularUniqueHarmonic::getK, &AngularUniqueHarmonic::setK)
    	;

      typedef class FixedTripleAngleListInteractionTemplate<AngularUniqueHarmonic>
        FixedTripleAngleListAngularUniqueHarmonic;
        
      class_ <FixedTripleAngleListAngularUniqueHarmonic, bases <Interaction> >
        ("interaction_FixedTripleAngleListAngularUniqueHarmonic",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleAngleList>,
                shared_ptr<AngularUniqueHarmonic> >())
        .def("setPotential", &FixedTripleAngleListAngularUniqueHarmonic::setPotential);
      ;
    }
  }
}
