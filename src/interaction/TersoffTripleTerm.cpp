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
#include "TersoffTripleTerm.hpp"

#include "VerletListTripleInteractionTemplate.hpp"
#include "FixedTripleListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    typedef class VerletListTripleInteractionTemplate <TersoffTripleTerm>
        VerletListTersoffTripleTerm;
    typedef class FixedTripleListInteractionTemplate <TersoffTripleTerm>
        FixedListTersoffTripleTerm;
    
    void 
    TersoffTripleTerm::registerPython() {
      using namespace espressopp::python;

      class_< TersoffTripleTerm, bases< AngularPotential > >(
            "interaction_TersoffTripleTerm",
            init< real, real, real, real, real, real, real, real, real, real, real, real, real, real >()
      )
      .add_property("B", &TersoffTripleTerm::getB,
              &TersoffTripleTerm::setB)
      .add_property("lambda2", &TersoffTripleTerm::getLambda2,
              &TersoffTripleTerm::setLambda2)
      .add_property("R", &TersoffTripleTerm::getR,
              &TersoffTripleTerm::setR)
      .add_property("D", &TersoffTripleTerm::getD,
              &TersoffTripleTerm::setD)
      .add_property("n", &TersoffTripleTerm::getN,
              &TersoffTripleTerm::setN)
      .add_property("beta", &TersoffTripleTerm::getBeta,
              &TersoffTripleTerm::setBeta)
      .add_property("m", &TersoffTripleTerm::getM,
              &TersoffTripleTerm::setM)
      .add_property("lambda3", &TersoffTripleTerm::getLambda3,
              &TersoffTripleTerm::setLambda3)
      .add_property("gamma", &TersoffTripleTerm::getGamma,
              &TersoffTripleTerm::setGamma)
      .add_property("c", &TersoffTripleTerm::getC,
              &TersoffTripleTerm::setC)
      .add_property("d", &TersoffTripleTerm::getd,
              &TersoffTripleTerm::setd)
      .add_property("theta0", &TersoffTripleTerm::getTheta0,
              &TersoffTripleTerm::setTheta0)
      .add_property("cutoff1", &TersoffTripleTerm::getCutoff1,
              &TersoffTripleTerm::setCutoff1)
      .add_property("cutoff2", &TersoffTripleTerm::getCutoff2,
              &TersoffTripleTerm::setCutoff2)
      ;

      
      class_< VerletListTersoffTripleTerm, bases< Interaction > >(
        "interaction_VerletListTersoffTripleTerm",
        init< shared_ptr<System>, 
              shared_ptr<VerletListTriple> >()
      )
      .def("getVerletListTriple", &VerletListTersoffTripleTerm::getVerletListTriple)
      .def("setPotential", &VerletListTersoffTripleTerm::setPotential,
              return_value_policy< reference_existing_object >())
      .def("getPotential", &VerletListTersoffTripleTerm::getPotential,
              return_value_policy< reference_existing_object >())
      ;

      class_< FixedListTersoffTripleTerm, bases< Interaction > >(
          "interaction_FixedTripleListTersoffTripleTerm",
          init<shared_ptr<System>, shared_ptr<FixedTripleList>, shared_ptr<TersoffTripleTerm> >()
      )
      .def("setPotential", &FixedListTersoffTripleTerm::setPotential)
      .def("getFixedTripleList", &FixedListTersoffTripleTerm::getFixedTripleList)
      ;
    }
  }
}
