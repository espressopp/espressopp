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
#include "GravityTruncated.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <GravityTruncated> VerletListGravityTruncated;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void GravityTruncated::registerPython() {
      using namespace espressopp::python;

      class_< GravityTruncated, bases< Potential > >
        ("interaction_GravityTruncated", init< >())
        .def(init< real, real >())
        .add_property("prefactor", &GravityTruncated::getPrefactor, &GravityTruncated::setPrefactor)
      ;

      class_< VerletListGravityTruncated, bases< Interaction > >
        ("interaction_VerletListGravityTruncated", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListGravityTruncated::getVerletList)
        .def("setPotential", &VerletListGravityTruncated::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListGravityTruncated::getPotential, return_value_policy< reference_existing_object >())
      ;
    }
  }
}
