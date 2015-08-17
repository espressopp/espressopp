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
#include "LennardJonesGenericCapped.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesGenericCapped>
        VerletListLennardJonesGenericCapped;
    typedef class VerletListAdressInteractionTemplate <LennardJonesGenericCapped, Tabulated>
        VerletListAdressLennardJonesGenericCapped;
    typedef class VerletListHadressInteractionTemplate <LennardJonesGenericCapped, Tabulated>
        VerletListHadressLennardJonesGenericCapped;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesGenericCapped> 
        CellListLennardJonesGenericCapped;
    typedef class FixedPairListInteractionTemplate <LennardJonesGenericCapped> 
        FixedPairListLennardJonesGenericCapped;
    LOG4ESPP_LOGGER(LennardJonesGenericCapped::theLogger, "LennardJonesGenericCapped");
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesGenericCapped::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesGenericCapped, bases< Potential > >
    	("interaction_LennardJonesGenericCapped", init< real, real, int, int, real, real, real >())
	    .def(init< real, real, int, int, real, real >())
    	.add_property("sigma", &LennardJonesGenericCapped::getSigma, &LennardJonesGenericCapped::setSigma)
    	.add_property("epsilon", &LennardJonesGenericCapped::getEpsilon, &LennardJonesGenericCapped::setEpsilon)
    	.add_property("a", &LennardJonesGenericCapped::getA, &LennardJonesGenericCapped::setA)
    	.add_property("b", &LennardJonesGenericCapped::getB, &LennardJonesGenericCapped::setB)
    	.add_property("caprad", &LennardJonesGenericCapped::getCaprad, &LennardJonesGenericCapped::setCaprad)
        .def_pickle(LennardJonesGenericCapped_pickle())

      ;

      class_< VerletListLennardJonesGenericCapped, bases< Interaction > > 
        ("interaction_VerletListLennardJonesGenericCapped", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesGenericCapped::getVerletList)
        .def("setPotential", &VerletListLennardJonesGenericCapped::setPotential)
        .def("getPotential", &VerletListLennardJonesGenericCapped::getPotentialPtr)
      ;

      class_< VerletListAdressLennardJonesGenericCapped, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesGenericCapped",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesGenericCapped::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesGenericCapped::setPotentialCG);
      ;
      
      class_< VerletListHadressLennardJonesGenericCapped, bases< Interaction > >
        ("interaction_VerletListHadressLennardJonesGenericCapped",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJonesGenericCapped::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJonesGenericCapped::setPotentialCG);
      ;
      
      class_< CellListLennardJonesGenericCapped, bases< Interaction > > 
        ("interaction_CellListLennardJonesGenericCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesGenericCapped::setPotential);
	  ;

      class_< FixedPairListLennardJonesGenericCapped, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesGenericCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesGenericCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesGenericCapped> >())
          .def("setPotential", &FixedPairListLennardJonesGenericCapped::setPotential)
          .def("getPotential", &FixedPairListLennardJonesGenericCapped::getPotential)
          .def("setFixedPairList", &FixedPairListLennardJonesGenericCapped::setFixedPairList)
          .def("getFixedPairList", &FixedPairListLennardJonesGenericCapped::getFixedPairList)
      ;
    }
    
  }
}
