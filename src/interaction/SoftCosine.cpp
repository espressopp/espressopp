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
#include "SoftCosine.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< SoftCosine > 
    VerletListSoftCosine;
    typedef class CellListAllPairsInteractionTemplate< SoftCosine > 
    CellListSoftCosine;
    typedef class FixedPairListInteractionTemplate< SoftCosine > 
    FixedPairListSoftCosine;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    SoftCosine::registerPython() {
      using namespace espressopp::python;

      class_< SoftCosine, bases< Potential > >
    	("interaction_SoftCosine", init< real, real >())
	.def(init< real, real, real>())
    	.add_property("A", &SoftCosine::getA, &SoftCosine::setA)
    	.def_pickle(SoftCosine_pickle())
    	;

      class_< VerletListSoftCosine, bases< Interaction > > 
        ("interaction_VerletListSoftCosine", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListSoftCosine::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListSoftCosine::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_< CellListSoftCosine, bases< Interaction > > 
        ("interaction_CellListSoftCosine", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListSoftCosine::setPotential);
	;

      class_< FixedPairListSoftCosine, bases< Interaction > >
        ("interaction_FixedPairListSoftCosine",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, 
                shared_ptr<SoftCosine> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<SoftCosine> >())
        .def("setPotential", &FixedPairListSoftCosine::setPotential);
        ;
    }
    
  }
}
