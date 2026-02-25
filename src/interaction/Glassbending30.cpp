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
#include "Glassbending30.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Glassbending30::registerPython() {
      using namespace espressopp::python;

      class_<Glassbending30, bases<AngularPotential> >("interaction_Glassbending30", init<real, real, real>())
        .add_property("Ba", &Glassbending30::getBa, &Glassbending30::setBa)
        .add_property("Bb", &Glassbending30::getBb, &Glassbending30::setBb)
        .add_property("Bbthreshold", &Glassbending30::getBbthreshold, &Glassbending30::setBbthreshold);


      typedef class FixedTripleListInteractionTemplate< Glassbending30 >
        FixedTripleListGlassbending30;
      class_< FixedTripleListGlassbending30, bases< Interaction > >
        ("interaction_FixedTripleListGlassbending30",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleList>,
                shared_ptr<Glassbending30> >())
        .def("setPotential", &FixedTripleListGlassbending30::setPotential)
        .def("getFixedTripleList", &FixedTripleListGlassbending30::getFixedTripleList);
      ;
    }
  }
}
