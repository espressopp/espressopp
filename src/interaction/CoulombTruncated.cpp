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
#include "CoulombTruncated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< CoulombTruncated > 
    VerletListCoulombTruncated;
    typedef class CellListAllPairsInteractionTemplate< CoulombTruncated > 
    CellListCoulombTruncated;
    typedef class FixedPairListInteractionTemplate< CoulombTruncated >
    FixedPairListCoulombTruncated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    CoulombTruncated::registerPython() {
      using namespace espressopp::python;

      class_< CoulombTruncated, bases< Potential > >
    	("interaction_CoulombTruncated", init< real, real >())
	    .def(init< real, real, real>())
    	.add_property("qq", &CoulombTruncated::getQQ, &CoulombTruncated::setQQ)
    	.def_pickle(CoulombTruncated_pickle())
    	;

      class_< VerletListCoulombTruncated, bases< Interaction > > 
        ("interaction_VerletListCoulombTruncated", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListCoulombTruncated::setPotential)
        .def("getPotential", &VerletListCoulombTruncated::getPotentialPtr)
        ;

      class_< CellListCoulombTruncated, bases< Interaction > > 
        ("interaction_CellListCoulombTruncated", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListCoulombTruncated::setPotential)
	;

      class_< FixedPairListCoulombTruncated, bases< Interaction > >
        ("interaction_FixedPairListCoulombTruncated",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<CoulombTruncated> >())
        .def("setPotential", &FixedPairListCoulombTruncated::setPotential)
        ;
    }
    
  }
}
