/*
  Copyright (C) 2017
      Max Planck Institute for Polymer Research
  
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
#include "ConstrainRG.hpp"
#include "FixedLocalTupleRgListInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {
	ConstrainRG::ConstrainRG(real _k_rg): k_rg(_k_rg){};
	
	ConstrainRG::~ConstrainRG(){};
	//////////////////////////////////////////////////
	// REGISTRATION WITH PYTHON
	//////////////////////////////////////////////////
	void ConstrainRG::registerPython() {
	    using namespace espressopp::python;
	    
	    class_< ConstrainRG, bases< Potential > >
		("interaction_ConstrainRG", init< real >() )
		.add_property("k_rg",
			      &ConstrainRG::getK_rg, 
			      &ConstrainRG::setK_rg)
		;
	    
	    typedef class FixedLocalTupleRgListInteractionTemplate<ConstrainRG>
		FixedLocalTupleListConstrainRG;
	    
	    class_< FixedLocalTupleListConstrainRG, bases< Interaction > >
		("interaction_FixedLocalTupleListConstrainRG",
		 init< shared_ptr<System>, shared_ptr<FixedLocalTupleList>, shared_ptr<ConstrainRG> >())
		.def("getPotential", &FixedLocalTupleListConstrainRG::getPotential)
		.def("setRG", &FixedLocalTupleListConstrainRG::setRG)
		;
	    
	}
	
    }
}
