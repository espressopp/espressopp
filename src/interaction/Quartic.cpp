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
#include "Quartic.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Quartic::registerPython() {
      using namespace espressopp::python;

      class_< Quartic, bases< Potential > >
    	("interaction_Quartic", init< real, real, real >())
	.def(init< real, real, real, real >())
	.add_property("K", &Quartic::getK, &Quartic::setK)
	.add_property("r0", &Quartic::getR0, &Quartic::setR0)
    	;

      typedef class FixedPairListInteractionTemplate< Quartic >
        FixedPairListQuartic;
      class_< FixedPairListQuartic, bases< Interaction > >
        ("interaction_FixedPairListQuartic",
           init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Quartic> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Quartic> >())
        .def("setPotential", &FixedPairListQuartic::setPotential)
        .def("getPotential", &FixedPairListQuartic::getPotential)
        .def("setFixedPairList", &FixedPairListQuartic::setFixedPairList)
        .def("getFixedPairList", &FixedPairListQuartic::getFixedPairList);
     ;
    }

  }
}
