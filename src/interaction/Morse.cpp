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
#include "Morse.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< Morse >
    VerletListMorse;
    typedef class VerletListAdressInteractionTemplate< Morse, Tabulated >
    VerletListAdressMorse;
    typedef class VerletListHadressInteractionTemplate< Morse, Tabulated >
    VerletListHadressMorse;
    typedef class CellListAllPairsInteractionTemplate< Morse >
    CellListMorse;
    typedef class FixedPairListInteractionTemplate< Morse >
    FixedPairListMorse;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Morse::registerPython() {
      using namespace espressopp::python;

      class_< Morse, bases< Potential > >
    	("interaction_Morse", init< real, real, real, real >())
	.def(init< real, real, real, real, real >())
    	.add_property("epsilon", &Morse::getEpsilon, &Morse::setEpsilon)
    	.add_property("alpha", &Morse::getAlpha, &Morse::setAlpha)
    	.add_property("rMin", &Morse::getRMin, &Morse::setRMin)
    	.def_pickle(Morse_pickle())
    	;

      class_< VerletListMorse, bases< Interaction > > 
        ("interaction_VerletListMorse", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListMorse::setPotential)
        .def("getPotential", &VerletListMorse::getPotentialPtr)
      ;

      class_< VerletListAdressMorse, bases< Interaction > >
        ("interaction_VerletListAdressMorse",
                init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressMorse::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressMorse::setPotentialCG);
        ;
        
      class_< VerletListHadressMorse, bases< Interaction > >
        ("interaction_VerletListHadressMorse",
                init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListHadressMorse::setPotentialAT)
        .def("setPotentialCG", &VerletListHadressMorse::setPotentialCG);
        ;
        
      class_< CellListMorse, bases< Interaction > > 
        ("interaction_CellListMorse", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListMorse::setPotential);
	;

      class_< FixedPairListMorse, bases< Interaction > >
        ("interaction_FixedPairListMorse",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Morse> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Morse> >())
        .def("setPotential", &FixedPairListMorse::setPotential);
        ;
    }
    
  }
}
