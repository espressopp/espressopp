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
#include "StillingerWeberPairTerm.hpp"
#include "Tabulated.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <StillingerWeberPairTerm>
        VerletListStillingerWeberPairTerm;
    typedef class VerletListAdressInteractionTemplate <StillingerWeberPairTerm, Tabulated>
        VerletListAdressStillingerWeberPairTerm;
    typedef class VerletListHadressInteractionTemplate <StillingerWeberPairTerm, Tabulated>
        VerletListHadressStillingerWeberPairTerm;
    typedef class CellListAllPairsInteractionTemplate <StillingerWeberPairTerm> 
        CellListStillingerWeberPairTerm;
    typedef class FixedPairListInteractionTemplate <StillingerWeberPairTerm> 
        FixedPairListStillingerWeberPairTerm;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    StillingerWeberPairTerm::registerPython() {
      using namespace espressopp::python;

      class_< StillingerWeberPairTerm, bases< Potential > >
    	("interaction_StillingerWeberPairTerm", init< real, real, real, real, real, real, real >())
	    //.def(init< real, real, real, real, real, real, real, real>())
        .add_property("A", &StillingerWeberPairTerm::getA, &StillingerWeberPairTerm::setA)
        .add_property("B", &StillingerWeberPairTerm::getB, &StillingerWeberPairTerm::setB)
        .add_property("p", &StillingerWeberPairTerm::getP, &StillingerWeberPairTerm::setP)
        .add_property("q", &StillingerWeberPairTerm::getQ, &StillingerWeberPairTerm::setQ)
    	.add_property("sigma", &StillingerWeberPairTerm::getSigma, &StillingerWeberPairTerm::setSigma)
    	.add_property("epsilon", &StillingerWeberPairTerm::getEpsilon, &StillingerWeberPairTerm::setEpsilon)
      ;

      class_< VerletListStillingerWeberPairTerm, bases< Interaction > > 
        ("interaction_VerletListStillingerWeberPairTerm", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListStillingerWeberPairTerm::getVerletList)
        .def("setPotential", &VerletListStillingerWeberPairTerm::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListStillingerWeberPairTerm::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressStillingerWeberPairTerm, bases< Interaction > >
        ("interaction_VerletListAdressStillingerWeberPairTerm",
            init< shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressStillingerWeberPairTerm::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressStillingerWeberPairTerm::setPotentialCG);
      ;

      class_< VerletListHadressStillingerWeberPairTerm, bases< Interaction > >
        ("interaction_VerletListHadressStillingerWeberPairTerm",
            init< shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressStillingerWeberPairTerm::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressStillingerWeberPairTerm::setPotentialCG);
      ;
      
      class_< CellListStillingerWeberPairTerm, bases< Interaction > > 
        ("interaction_CellListStillingerWeberPairTerm", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListStillingerWeberPairTerm::setPotential);
	  ;

      class_< FixedPairListStillingerWeberPairTerm, bases< Interaction > >
        ("interaction_FixedPairListStillingerWeberPairTerm",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<StillingerWeberPairTerm> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<StillingerWeberPairTerm> >())
          .def("setPotential", &FixedPairListStillingerWeberPairTerm::setPotential);
      ;
    }
    
  }
}
