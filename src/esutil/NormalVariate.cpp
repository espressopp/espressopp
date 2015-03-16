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
#include "NormalVariate.hpp"
#include "esutil/RNG.hpp"

using namespace boost;

namespace espressopp {
  namespace esutil {
    NormalVariate::NormalVariate(shared_ptr< RNG > _rng,
				 const real mean, 
				 const real sigma) 
      : Super(*(_rng->getBoostRNG()), DistType(mean, sigma)), rng(_rng)
    {}
    
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    NormalVariate::registerPython() {
      using namespace espressopp::python;

      real (NormalVariate::*pyCall)() 
      	= &NormalVariate::operator();

      class_< NormalVariate >("esutil_NormalVariate",
			      init< shared_ptr< RNG > >())
      	.def("__call__", pyCall)
      	;
    }
  }
}
