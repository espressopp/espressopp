/*
  Copyright (C) 2012,2013,2016
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
#include "LennardJonesSoftcoreTI.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesSoftcoreTI>
        VerletListLennardJonesSoftcoreTI;
    typedef class VerletListAdressInteractionTemplate <LennardJonesSoftcoreTI, Tabulated>
        VerletListAdressLennardJonesSoftcoreTI;
    typedef class VerletListHadressInteractionTemplate <LennardJonesSoftcoreTI, Tabulated>
        VerletListHadressLennardJonesSoftcoreTI;
    typedef class FixedPairListInteractionTemplate <LennardJonesSoftcoreTI> 
        FixedPairListLennardJonesSoftcoreTI;
    typedef class FixedPairListTypesInteractionTemplate <LennardJonesSoftcoreTI> 
        FixedPairListTypesLennardJonesSoftcoreTI;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesSoftcoreTI::registerPython() {
      using namespace espressopp::python;

      void (LennardJonesSoftcoreTI::*pyAddPid)(longint pid) = &LennardJonesSoftcoreTI::addPid;

      class_< LennardJonesSoftcoreTI, bases< Potential > >
    	("interaction_LennardJonesSoftcoreTI", init< real, real, real, real, real, real, real, real, bool>())
        .def("addPid",pyAddPid)
      ;

      //class_< VerletListLennardJonesSoftcoreTI, bases< Interaction > > 
      //  ("interaction_VerletListLennardJonesSoftcoreTI", init< shared_ptr<VerletList> >())
      //  .def("getVerletList", &VerletListLennardJonesSoftcoreTI::getVerletList)
      //  .def("setPotential", &VerletListLennardJonesSoftcoreTI::setPotential)
      //  .def("getPotential", &VerletListLennardJonesSoftcoreTI::getPotentialPtr)
      //;

      class_< VerletListAdressLennardJonesSoftcoreTI, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesSoftcoreTI",
           init< shared_ptr<VerletListAdress>,
                  shared_ptr<FixedTupleListAdress> >())
        .def("setPotentialAT", &VerletListAdressLennardJonesSoftcoreTI::setPotentialAT)
        .def("setPotentialCG", &VerletListAdressLennardJonesSoftcoreTI::setPotentialCG);
      ;
      
      //class_< VerletListHadressLennardJonesSoftcoreTI, bases< Interaction > >
      //  ("interaction_VerletListHadressLennardJonesSoftcoreTI",
      //     init< shared_ptr<VerletListAdress>,
      //            shared_ptr<FixedTupleListAdress> >())
      //  .def("setPotentialAT", &VerletListHadressLennardJonesSoftcoreTI::setPotentialAT)
      //  .def("setPotentialCG", &VerletListHadressLennardJonesSoftcoreTI::setPotentialCG);
      //;
      
    }
    
  }
}
