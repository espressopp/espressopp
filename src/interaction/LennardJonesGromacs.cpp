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
#include "LennardJonesGromacs.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    typedef class VerletListInteractionTemplate< LennardJonesGromacs >
    VerletListLennardJonesGromacs;
    typedef class CellListAllPairsInteractionTemplate< LennardJonesGromacs >
    CellListLennardJonesGromacs;
    typedef class FixedPairListInteractionTemplate< LennardJonesGromacs >
    FixedPairListLennardJonesGromacs;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesGromacs::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesGromacs, bases< Potential > >
    	("interaction_LennardJonesGromacs", init< real, real, real, real >())
	.def(init< real, real, real, real, real >())
    	.add_property("epsilon", &LennardJonesGromacs::getEpsilon, &LennardJonesGromacs::setEpsilon)
    	.add_property("sigma", &LennardJonesGromacs::getSigma, &LennardJonesGromacs::setSigma)
    	.add_property("r1", &LennardJonesGromacs::getR1, &LennardJonesGromacs::setR1)
    	.def_pickle(LennardJonesGromacs_pickle())
    	;

      class_< VerletListLennardJonesGromacs, bases< Interaction > > 
        ("interaction_VerletListLennardJonesGromacs", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesGromacs::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesGromacs::getPotential, return_value_policy< reference_existing_object >())
        ;

      class_< CellListLennardJonesGromacs, bases< Interaction > >
        ("interaction_CellListLennardJonesGromacs", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesGromacs::setPotential);
	;

      class_< FixedPairListLennardJonesGromacs, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesGromacs",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesGromacs> >())
        .def("setPotential", &FixedPairListLennardJonesGromacs::setPotential);
        ;
    }
  }
}
