/*
  Copyright (C) 2012,2013,2015
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
#include "CoulombTruncatedUniqueCharge.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< CoulombTruncatedUniqueCharge > 
    VerletListCoulombTruncatedUniqueCharge;
    typedef class CellListAllPairsInteractionTemplate< CoulombTruncatedUniqueCharge > 
    CellListCoulombTruncatedUniqueCharge;
    typedef class FixedPairListInteractionTemplate< CoulombTruncatedUniqueCharge >
    FixedPairListCoulombTruncatedUniqueCharge;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    CoulombTruncatedUniqueCharge::registerPython() {
      using namespace espressopp::python;

      class_< CoulombTruncatedUniqueCharge, bases< Potential > >
    	("interaction_CoulombTruncatedUniqueCharge", init< real, real >())
	    .def(init< real, real, real>())
    	.add_property("qq", &CoulombTruncatedUniqueCharge::getQQ, &CoulombTruncatedUniqueCharge::setQQ)
    	.def_pickle(CoulombTruncatedUniqueCharge_pickle())
    	;

      class_< VerletListCoulombTruncatedUniqueCharge, bases< Interaction > > 
        ("interaction_VerletListCoulombTruncatedUniqueCharge", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListCoulombTruncatedUniqueCharge::setPotential)
        .def("getPotential", &VerletListCoulombTruncatedUniqueCharge::getPotentialPtr)
        ;

      class_< CellListCoulombTruncatedUniqueCharge, bases< Interaction > > 
        ("interaction_CellListCoulombTruncatedUniqueCharge", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListCoulombTruncatedUniqueCharge::setPotential)
	;

      class_< FixedPairListCoulombTruncatedUniqueCharge, bases< Interaction > >
        ("interaction_FixedPairListCoulombTruncatedUniqueCharge",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<CoulombTruncatedUniqueCharge> >())
        .def("setPotential", &FixedPairListCoulombTruncatedUniqueCharge::setPotential)
        ;
    }
    
  }
}
