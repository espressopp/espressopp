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
#include "StillingerWeberPairTermCapped.hpp"
#include "Tabulated.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <StillingerWeberPairTermCapped>
        VerletListStillingerWeberPairTermCapped;
    typedef class VerletListAdressInteractionTemplate <StillingerWeberPairTermCapped, Tabulated>
        VerletListAdressStillingerWeberPairTermCapped;
    typedef class VerletListHadressInteractionTemplate <StillingerWeberPairTermCapped, Tabulated>
        VerletListHadressStillingerWeberPairTermCapped;
    typedef class CellListAllPairsInteractionTemplate <StillingerWeberPairTermCapped> 
        CellListStillingerWeberPairTermCapped;
    typedef class FixedPairListInteractionTemplate <StillingerWeberPairTermCapped> 
        FixedPairListStillingerWeberPairTermCapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    StillingerWeberPairTermCapped::registerPython() {
      using namespace espressopp::python;

      class_< StillingerWeberPairTermCapped, bases< Potential > >
    	("interaction_StillingerWeberPairTermCapped",
              init< real, real, real, real, real, real, real, real >())
        .add_property("A", &StillingerWeberPairTermCapped::getA, &StillingerWeberPairTermCapped::setA)
        .add_property("B", &StillingerWeberPairTermCapped::getB, &StillingerWeberPairTermCapped::setB)
        .add_property("p", &StillingerWeberPairTermCapped::getP, &StillingerWeberPairTermCapped::setP)
        .add_property("q", &StillingerWeberPairTermCapped::getQ, &StillingerWeberPairTermCapped::setQ)
    	.add_property("sigma", &StillingerWeberPairTermCapped::getSigma, &StillingerWeberPairTermCapped::setSigma)
    	.add_property("epsilon", &StillingerWeberPairTermCapped::getEpsilon, &StillingerWeberPairTermCapped::setEpsilon)
    	.add_property("caprad", &StillingerWeberPairTermCapped::getCaprad, &StillingerWeberPairTermCapped::setCaprad)
        .def("getCaprad", &StillingerWeberPairTermCapped::getCaprad)
      ;

      class_< VerletListStillingerWeberPairTermCapped, bases< Interaction > > 
        ("interaction_VerletListStillingerWeberPairTermCapped", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListStillingerWeberPairTermCapped::getVerletList)
        .def("setPotential", &VerletListStillingerWeberPairTermCapped::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListStillingerWeberPairTermCapped::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressStillingerWeberPairTermCapped, bases< Interaction > >
        ("interaction_VerletListAdressStillingerWeberPairTermCapped",
            init< shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressStillingerWeberPairTermCapped::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressStillingerWeberPairTermCapped::setPotentialCG);
      ;
      
      class_< VerletListHadressStillingerWeberPairTermCapped, bases< Interaction > >
        ("interaction_VerletListHadressStillingerWeberPairTermCapped",
            init< shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressStillingerWeberPairTermCapped::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressStillingerWeberPairTermCapped::setPotentialCG);
      ;

      class_< CellListStillingerWeberPairTermCapped, bases< Interaction > > 
        ("interaction_CellListStillingerWeberPairTermCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListStillingerWeberPairTermCapped::setPotential);
	  ;

      class_< FixedPairListStillingerWeberPairTermCapped, bases< Interaction > >
        ("interaction_FixedPairListStillingerWeberPairTermCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<StillingerWeberPairTermCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<StillingerWeberPairTermCapped> >())
          .def("setPotential", &FixedPairListStillingerWeberPairTermCapped::setPotential);
      ;
    }
    
  }
}
