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
#include "LennardJonesExpand.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< LennardJonesExpand >
    VerletListLennardJonesExpand;
    typedef class CellListAllPairsInteractionTemplate< LennardJonesExpand >
    CellListLennardJonesExpand;
    typedef class FixedPairListInteractionTemplate< LennardJonesExpand >
    FixedPairListLennardJonesExpand;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesExpand::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesExpand, bases< Potential > >
    	("interaction_LennardJonesExpand", init< real, real, real, real >())
	.def(init< real, real, real, real, real >())
    	.add_property("epsilon", &LennardJonesExpand::getEpsilon, &LennardJonesExpand::setEpsilon)
    	.add_property("sigma", &LennardJonesExpand::getSigma, &LennardJonesExpand::setSigma)
    	.add_property("delta", &LennardJonesExpand::getDelta, &LennardJonesExpand::setDelta)
    	.def_pickle(LennardJonesExpand_pickle())
    	;

      class_< VerletListLennardJonesExpand, bases< Interaction > > 
        ("interaction_VerletListLennardJonesExpand", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesExpand::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesExpand::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_< CellListLennardJonesExpand, bases< Interaction > >
        ("interaction_CellListLennardJonesExpand", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesExpand::setPotential);
	;

      class_< FixedPairListLennardJonesExpand, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesExpand",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesExpand> >())
        .def("setPotential", &FixedPairListLennardJonesExpand::setPotential);
        ;
    }
  }
}
