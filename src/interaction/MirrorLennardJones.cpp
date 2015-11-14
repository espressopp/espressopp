/*
  Copyright (C) 2014
      Pierre de Buyl
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
#include "MirrorLennardJones.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
      typedef class FixedPairListInteractionTemplate< MirrorLennardJones >
      FixedPairListMirrorLennardJones;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    MirrorLennardJones::registerPython() {
      using namespace espressopp::python;

      class_< MirrorLennardJones, bases< Potential > >
    	("interaction_MirrorLennardJones", init< real, real >())
	.def(init< real, real >())
	.add_property("epsilon", &MirrorLennardJones::getEpsilon, &MirrorLennardJones::setEpsilon)
	.add_property("sigma", &MirrorLennardJones::getSigma, &MirrorLennardJones::setSigma)
    	;

      class_< FixedPairListMirrorLennardJones, bases< Interaction > >
      ("interaction_FixedPairListMirrorLennardJones",
        init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<MirrorLennardJones> >())
       .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<MirrorLennardJones> >())
       .def("setPotential", &FixedPairListMirrorLennardJones::setPotential)
       .def("getPotential", &FixedPairListMirrorLennardJones::getPotential)
       .def("setFixedPairList", &FixedPairListMirrorLennardJones::setFixedPairList)
       .def("getFixedPairList", &FixedPairListMirrorLennardJones::getFixedPairList)
       ;
    }

  }
}
