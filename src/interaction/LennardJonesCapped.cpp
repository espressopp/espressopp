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
#include "LennardJonesCapped.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesCapped>
        VerletListLennardJonesCapped;
    typedef class VerletListAdressInteractionTemplate <LennardJonesCapped, Tabulated>
        VerletListAdressLennardJonesCapped;
    typedef class VerletListHadressInteractionTemplate <LennardJonesCapped, Tabulated>
        VerletListHadressLennardJonesCapped;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesCapped>
        CellListLennardJonesCapped;
    typedef class FixedPairListInteractionTemplate <LennardJonesCapped>
        FixedPairListLennardJonesCapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesCapped::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesCapped, bases< Potential > >
    	("interaction_LennardJonesCapped", init< real, real, real, real >())
	    .def(init< real, real, real, real, real>())
    	.add_property("sigma", &LennardJonesCapped::getSigma, &LennardJonesCapped::setSigma)
    	.add_property("epsilon", &LennardJonesCapped::getEpsilon, &LennardJonesCapped::setEpsilon)
    	.add_property("caprad", &LennardJonesCapped::getCaprad, &LennardJonesCapped::setCaprad)
    	.def_pickle(LennardJonesCapped_pickle())
      ;

      class_< VerletListLennardJonesCapped, bases< Interaction > >
        ("interaction_VerletListLennardJonesCapped", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesCapped::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesCapped::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressLennardJonesCapped, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesCapped",
         init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
         .def("setPotentialAT", &VerletListAdressLennardJonesCapped::setPotentialAT)
         .def("setPotentialCG", &VerletListAdressLennardJonesCapped::setPotentialCG)
         .def("getPotentialAT", &VerletListAdressLennardJonesCapped::getPotentialAT,
                 return_value_policy< reference_existing_object >())
         .def("getPotentialCG", &VerletListAdressLennardJonesCapped::getPotentialCG,
                       return_value_policy< reference_existing_object >());
      ;

      class_< VerletListHadressLennardJonesCapped, bases< Interaction > >
        ("interaction_VerletListHadressLennardJonesCapped",
         init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
         .def("setPotentialAT", &VerletListHadressLennardJonesCapped::setPotentialAT)
         .def("setPotentialCG", &VerletListHadressLennardJonesCapped::setPotentialCG)
         .def("getPotentialAT", &VerletListHadressLennardJonesCapped::getPotentialAT,
                 return_value_policy< reference_existing_object >())
         .def("getPotentialCG", &VerletListHadressLennardJonesCapped::getPotentialCG,
                       return_value_policy< reference_existing_object >());
      ;
      
      class_< CellListLennardJonesCapped, bases< Interaction > >
        ("interaction_CellListLennardJonesCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesCapped::setPotential)
        .def("getPotential", &CellListLennardJonesCapped::getPotential, return_value_policy< reference_existing_object >());
	  ;

      class_< FixedPairListLennardJonesCapped, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesCapped> >())
          .def("setPotential", &FixedPairListLennardJonesCapped::setPotential);
      ;
    }
    
  }
}
