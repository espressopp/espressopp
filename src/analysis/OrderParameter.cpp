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
#include "OrderParameter.hpp"

#include "boost/math/special_functions/spherical_harmonic.hpp"

namespace espressopp {
  namespace analysis {


    dcomplex OrderParameter::SphHarm(int l, int m, Real3D r) {
      real d = r.abs(); // distance between two particles
      real theta, phi;

      // defining theta and phi
      theta = acos( r[2] / d );   // in radians
      
      // phi is defined as http://en.wikipedia.org/wiki/Atan2
      // problem is x = y = 0, we will define it like 0
      if(r[0] > 0.0){
        phi = atan( r[1]/r[0] );
      }
      else if( r[0] < 0.0 && r[1] >= 0.0 ){
        phi = atan( r[1]/r[0] ) + M_PIl;
      }
      else if( r[0] < 0.0 && r[1] < 0.0 ){
        phi = atan( r[1]/r[0] ) - M_PIl;
      }
      else if( r[0] == 0.0 && r[1] > 0.0 ){
        phi = M_PIl;
      }
      else if( r[0] == 0.0 && r[1] < 0.0 ){
        phi = -M_PIl;
      }
      else{
        // x = 0; y = 0;
        phi = 0.0;
      }
      
      return boost::math::spherical_harmonic(l, m, theta, phi);
      //throw std::runtime_error("OrderParameter::SphHarm is broken (used for coupled cluster analysis only");
      //return dcomplex(0);
    }
    
    void OrderParameter::registerPython() {
      using namespace espressopp::python;
      class_<OrderParameter, bases< AnalysisBase > >
        ("analysis_OrderParameter", init< shared_ptr<System>, real, int, bool, bool, real, real>())
        .add_property("l", 
              &OrderParameter::getAngularMomentum,
              &OrderParameter::setAngularMomentum)
        .add_property("cutoff", 
              &OrderParameter::getCutoff,
              &OrderParameter::setCutoff)
      ;
    }
  }
}
