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
#include "CoulombRSpace.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"

// currently just Verlet list

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <CoulombRSpace> VerletListCoulombRSpace;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void CoulombRSpace::registerPython() {
      using namespace espressopp::python;

      class_< CoulombRSpace, bases< Potential > >
        ("interaction_CoulombRSpace", init< >())
        .def(init< real, real, real >())
        .add_property("alpha", &CoulombRSpace::getAlpha, &CoulombRSpace::setAlpha)
        .add_property("prefactor", &CoulombRSpace::getPrefactor, &CoulombRSpace::setPrefactor)
      ;

      class_< VerletListCoulombRSpace, bases< Interaction > >
        ("interaction_VerletListCoulombRSpace", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListCoulombRSpace::getVerletList)
        .def("setPotential", &VerletListCoulombRSpace::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListCoulombRSpace::getPotential, return_value_policy< reference_existing_object >())
      ;
    }
    
  }
}
