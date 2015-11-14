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
#include "LennardJonesGeneric.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesGeneric>
        VerletListLennardJonesGeneric;
    typedef class VerletListAdressInteractionTemplate <LennardJonesGeneric, Tabulated>
        VerletListAdressLennardJonesGeneric;
    typedef class VerletListAdressInteractionTemplate <LennardJonesGeneric, LennardJonesGeneric>
        VerletListAdressLennardJonesGeneric2;
    typedef class VerletListHadressInteractionTemplate <LennardJonesGeneric, Tabulated>
        VerletListHadressLennardJonesGeneric;
    typedef class VerletListHadressInteractionTemplate <LennardJonesGeneric, LennardJonesGeneric>
        VerletListHadressLennardJonesGeneric2;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesGeneric> 
        CellListLennardJonesGeneric;
    typedef class FixedPairListInteractionTemplate <LennardJonesGeneric> 
        FixedPairListLennardJonesGeneric;
    LOG4ESPP_LOGGER(LennardJonesGeneric::theLogger, "LennardJonesGeneric");
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesGeneric::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesGeneric, bases< Potential > >
    	("interaction_LennardJonesGeneric", init< real, real, int, int, real, real >())
	    .def(init< real, real, int, int, real >())
    	.add_property("sigma", &LennardJonesGeneric::getSigma, &LennardJonesGeneric::setSigma)
    	.add_property("epsilon", &LennardJonesGeneric::getEpsilon, &LennardJonesGeneric::setEpsilon)
    	.add_property("a", &LennardJonesGeneric::getA, &LennardJonesGeneric::setA)
    	.add_property("b", &LennardJonesGeneric::getB, &LennardJonesGeneric::setB)
        .def_pickle(LennardJonesGeneric_pickle())

      ;

      class_< VerletListLennardJonesGeneric, bases< Interaction > > 
        ("interaction_VerletListLennardJonesGeneric", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesGeneric::getVerletList)
        .def("setPotential", &VerletListLennardJonesGeneric::setPotential)
        .def("getPotential", &VerletListLennardJonesGeneric::getPotentialPtr)
      ;

      class_< VerletListAdressLennardJonesGeneric, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesGeneric",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesGeneric::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesGeneric::setPotentialCG);
      ;
      
      class_< VerletListAdressLennardJonesGeneric2, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesGeneric2",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesGeneric2::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesGeneric2::setPotentialCG);
      ;

      class_< VerletListHadressLennardJonesGeneric, bases< Interaction > >
        ("interaction_VerletListHadressLennardJonesGeneric",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJonesGeneric::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJonesGeneric::setPotentialCG);
      ;
      
      class_< VerletListHadressLennardJonesGeneric2, bases< Interaction > >
        ("interaction_VerletListHadressLennardJonesGeneric2",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressLennardJonesGeneric2::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressLennardJonesGeneric2::setPotentialCG);
      ;
      
      class_< CellListLennardJonesGeneric, bases< Interaction > > 
        ("interaction_CellListLennardJonesGeneric", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesGeneric::setPotential);
	  ;

      class_< FixedPairListLennardJonesGeneric, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesGeneric",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesGeneric> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesGeneric> >())
          .def("setPotential", &FixedPairListLennardJonesGeneric::setPotential)
          .def("getPotential", &FixedPairListLennardJonesGeneric::getPotential)
          .def("setFixedPairList", &FixedPairListLennardJonesGeneric::setFixedPairList)
          .def("getFixedPairList", &FixedPairListLennardJonesGeneric::getFixedPairList)
      ;
    }
    
  }
}
