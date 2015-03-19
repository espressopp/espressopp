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
#include "AngularUniqueCosineSquared.hpp"
#include "FixedTripleAngleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularUniqueCosineSquared::registerPython() {
      using namespace espressopp::python;

      class_< AngularUniqueCosineSquared, bases< AngularUniquePotential > >
    	("interaction_AngularUniqueCosineSquared", init< real >())
        .add_property("K", &AngularUniqueCosineSquared::getK, &AngularUniqueCosineSquared::setK)
      ;

      typedef class FixedTripleAngleListInteractionTemplate< AngularUniqueCosineSquared >
        FixedTripleAngleListAngularUniqueCosineSquared;
        class_< FixedTripleAngleListAngularUniqueCosineSquared, bases< Interaction > >
          ("interaction_FixedTripleAngleListAngularUniqueCosineSquared",
                  init<shared_ptr<System>,
                      shared_ptr<FixedTripleAngleList>,
                      shared_ptr<AngularUniqueCosineSquared> >())
          .def("setPotential", &FixedTripleAngleListAngularUniqueCosineSquared::setPotential)
          .def("getFixedTripleList", &FixedTripleAngleListAngularUniqueCosineSquared::getFixedTripleList)
        ;
    }
  }
}
