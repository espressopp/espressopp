/*
  Copyright (C) 2012-2016
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
#include "FENE.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class FixedPairListInteractionTemplate< FENE > FixedPairListFENE;
    LOG4ESPP_LOGGER(FENE::theLogger, "FENE");
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    FENE::registerPython() {
      using namespace espressopp::python;
      class_< FENE, bases< Potential > >
      ("interaction_FENE", init< real, real, real, real >())
	  .def(init< real, real, real, real, real >())
	  .add_property("K", &FENE::getK, &FENE::setK)
	  .add_property("r0", &FENE::getR0, &FENE::setR0)
	  .add_property("rMax", &FENE::getRMax, &FENE::setRMax)
	  ;

      class_< FixedPairListFENE, bases< Interaction > >
      ("interaction_FixedPairListFENE", init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<FENE> >())
      .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<FENE> >())
      .def("setPotential", &FixedPairListFENE::setPotential)
      .def("getPotential", &FixedPairListFENE::getPotential)
      .def("setFixedPairList", &FixedPairListFENE::setFixedPairList)
      .def("getFixedPairList", &FixedPairListFENE::getFixedPairList)
      ;
    }

  }
}
