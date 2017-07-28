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
#include "ConstrainCOM.hpp"
#include "FixedLocalTupleComListInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {
	ConstrainCOM::ConstrainCOM(real _k_com): k_com(_k_com){};
	
	ConstrainCOM::~ConstrainCOM(){};
	//////////////////////////////////////////////////
	// REGISTRATION WITH PYTHON
	//////////////////////////////////////////////////
	void ConstrainCOM::registerPython() {
	    using namespace espressopp::python;
	    
	    class_< ConstrainCOM, bases< Potential > >
		("interaction_ConstrainCOM", init< real >() )
		.add_property("k_com",
			      &ConstrainCOM::getK_com, 
			      &ConstrainCOM::setK_com)
		;
	    
	    typedef class FixedLocalTupleComListInteractionTemplate<ConstrainCOM>
		FixedLocalTupleListConstrainCOM;
	    
	    class_< FixedLocalTupleListConstrainCOM, bases< Interaction > >
		("interaction_FixedLocalTupleListConstrainCOM",
		 init< shared_ptr<System>, shared_ptr<FixedLocalTupleList>, shared_ptr<ConstrainCOM> >())
		.def("getPotential", &FixedLocalTupleListConstrainCOM::getPotential)
		.def("setCom", &FixedLocalTupleListConstrainCOM::setCom)
		;
	    
	}
	
    }
}
