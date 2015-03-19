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
#include "HarmonicUnique.hpp"
#include "FixedPairDistListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    HarmonicUnique::registerPython() {
      using namespace espressopp::python;

      class_< HarmonicUnique, bases< PotentialUniqueDist > >
    	("interaction_HarmonicUnique", init< real >())
        .add_property("K", &HarmonicUnique::getK, &HarmonicUnique::setK)
      ;

      typedef class FixedPairDistListInteractionTemplate< HarmonicUnique >
        FixedPairDistListHarmonicUnique;
      class_< FixedPairDistListHarmonicUnique, bases< Interaction > >
        ("interaction_FixedPairDistListHarmonicUnique",
           init< shared_ptr<System>, shared_ptr<FixedPairDistList>, shared_ptr<HarmonicUnique> >())
        .def("setPotential", &FixedPairDistListHarmonicUnique::setPotential)
        .def("getPotential", &FixedPairDistListHarmonicUnique::getPotential)
        .def("setFixedPairList", &FixedPairDistListHarmonicUnique::setFixedPairList)
        .def("getFixedPairList", &FixedPairDistListHarmonicUnique::getFixedPairList);
      ;
    }

  }
}
