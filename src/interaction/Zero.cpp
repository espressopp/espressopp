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
#include "Zero.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
      
    typedef class VerletListInteractionTemplate <Zero>
        VerletListZero;
    typedef class VerletListAdressInteractionTemplate <Zero, Tabulated>
        VerletListAdressZero;
    typedef class VerletListHadressInteractionTemplate <Zero, Tabulated>
        VerletListHadressZero;
    typedef class CellListAllPairsInteractionTemplate <Zero>
        CellListZero;
    typedef class FixedPairListInteractionTemplate <Zero>
        FixedPairListZero;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Zero::registerPython() {
      using namespace espressopp::python;

      class_< Zero, bases< Potential > >
    	("interaction_Zero", init<>())
	    .def(init<>())
	    .def_pickle(Zero_pickle())
      ;

      class_< VerletListZero, bases< Interaction > >
        ("interaction_VerletListZero", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListZero::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListZero::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressZero, bases< Interaction > >
        ("interaction_VerletListAdressZero",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setFixedTupleList", &VerletListAdressZero::setFixedTupleList)
        .def("setPotentialAT", &VerletListAdressZero::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressZero::setPotentialCG);
      ;

      class_< VerletListHadressZero, bases< Interaction > >
        ("interaction_VerletListHadressZero",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setFixedTupleList", &VerletListHadressZero::setFixedTupleList)
        .def("setPotentialAT", &VerletListHadressZero::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressZero::setPotentialCG);
      ;
      
      class_< CellListZero, bases< Interaction > >
        ("interaction_CellListZero", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListZero::setPotential);
	  ;

      class_< FixedPairListZero, bases< Interaction > >
        ("interaction_FixedPairListZero",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Zero> >())
          .def("setPotential", &FixedPairListZero::setPotential);
      ;
    }
    
  }
}
