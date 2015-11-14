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
#include "Cosine.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Cosine::registerPython() {
      using namespace espressopp::python;

      class_< Cosine, bases< AngularPotential > >
    	("interaction_Cosine", init< real, real >())
	.add_property("K", &Cosine::getK, &Cosine::setK)
	.add_property("theta0", &Cosine::getTheta0, &Cosine::setTheta0)
    	;

      typedef class FixedTripleListInteractionTemplate< Cosine >
        FixedTripleListCosine;
      class_< FixedTripleListCosine, bases< Interaction > >
        ("interaction_FixedTripleListCosine",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleList>,
                shared_ptr<Cosine> >())
        .def("setPotential", &FixedTripleListCosine::setPotential)
        .def("getFixedTripleList", &FixedTripleListCosine::getFixedTripleList);
      ;
    }
  }
}
