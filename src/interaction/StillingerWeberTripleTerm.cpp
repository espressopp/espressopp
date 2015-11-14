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
#include "StillingerWeberTripleTerm.hpp"

#include "VerletListTripleInteractionTemplate.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    typedef class VerletListTripleInteractionTemplate <StillingerWeberTripleTerm>
        VerletListStillingerWeberTripleTerm;
    typedef class FixedTripleListInteractionTemplate <StillingerWeberTripleTerm>
        FixedListStillingerWeberTripleTerm;
    
    void 
    StillingerWeberTripleTerm::registerPython() {
      using namespace espressopp::python;

      class_< StillingerWeberTripleTerm, bases< AngularPotential > >(
            "interaction_StillingerWeberTripleTerm",
            init< real, real, real, real, real, real, real, real, real >()
      )
      .add_property("gamma1", &StillingerWeberTripleTerm::getGamma1,
              &StillingerWeberTripleTerm::setGamma1)
      .add_property("gamma2", &StillingerWeberTripleTerm::getGamma2,
              &StillingerWeberTripleTerm::setGamma2)
      .add_property("theta0", &StillingerWeberTripleTerm::getTheta0,
              &StillingerWeberTripleTerm::setTheta0)
      .add_property("lambda", &StillingerWeberTripleTerm::getLambda,
              &StillingerWeberTripleTerm::setLambda)
      .add_property("epsilon", &StillingerWeberTripleTerm::getEpsilon,
              &StillingerWeberTripleTerm::setEpsilon)
      .add_property("sigma1", &StillingerWeberTripleTerm::getSigma1,
              &StillingerWeberTripleTerm::setSigma1)
      .add_property("sigma2", &StillingerWeberTripleTerm::getSigma2,
              &StillingerWeberTripleTerm::setSigma2)
      .add_property("cutoff1", &StillingerWeberTripleTerm::getCutoff1,
              &StillingerWeberTripleTerm::setCutoff1)
      .add_property("cutoff2", &StillingerWeberTripleTerm::getCutoff2,
              &StillingerWeberTripleTerm::setCutoff2)
      ;

      
      class_< VerletListStillingerWeberTripleTerm, bases< Interaction > >(
        "interaction_VerletListStillingerWeberTripleTerm",
        init< shared_ptr<System>, 
              shared_ptr<VerletListTriple> >()
      )
      .def("getVerletListTriple", &VerletListStillingerWeberTripleTerm::getVerletListTriple)
      .def("setPotential", &VerletListStillingerWeberTripleTerm::setPotential,
              return_value_policy< reference_existing_object >())
      .def("getPotential", &VerletListStillingerWeberTripleTerm::getPotential,
              return_value_policy< reference_existing_object >())
      ;

      class_< FixedListStillingerWeberTripleTerm, bases< Interaction > >(
          "interaction_FixedTripleListStillingerWeberTripleTerm",
          init<shared_ptr<System>, shared_ptr<FixedTripleList>, shared_ptr<StillingerWeberTripleTerm> >()
      )
      .def("setPotential", &FixedListStillingerWeberTripleTerm::setPotential)
      .def("getFixedTripleList", &FixedListStillingerWeberTripleTerm::getFixedTripleList)
      ;
    }
  }
}
