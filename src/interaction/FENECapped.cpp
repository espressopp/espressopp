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
#include "FENECapped.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
      typedef class FixedPairListInteractionTemplate< FENECapped >
      FixedPairListFENECapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    FENECapped::registerPython() {
      using namespace espressopp::python;

      class_< FENECapped, bases< Potential > >
    	("interaction_FENECapped", init< real, real, real, real, real >())
	.def(init< real, real, real, real, real, real >())
	.add_property("K", &FENECapped::getK, &FENECapped::setK)
   .add_property("r_cap", &FENECapped::getRCap, &FENECapped::setRCap)
	.add_property("r0", &FENECapped::getR0, &FENECapped::setR0)
	.add_property("rMax", &FENECapped::getRMax, &FENECapped::setRMax)
    	;

      class_< FixedPairListFENECapped, bases< Interaction > >
      ("interaction_FixedPairListFENECapped",
        init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<FENECapped> >())
       .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<FENECapped> >())
       .def("setPotential", &FixedPairListFENECapped::setPotential)
       .def("getPotential", &FixedPairListFENECapped::getPotential)
       .def("setFixedPairList", &FixedPairListFENECapped::setFixedPairList)
       .def("getFixedPairList", &FixedPairListFENECapped::getFixedPairList)
       ;
    }

  }
}
