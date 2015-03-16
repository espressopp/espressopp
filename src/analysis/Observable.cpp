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

#include <python.hpp>
#include "Observable.hpp"
#include "boost/foreach.hpp"

namespace espressopp {
  namespace analysis {

    LOG4ESPP_LOGGER(Observable::logger, "Observable");

    python::list Observable::compute_real_vector_python() {
      python::list ret;
      compute_real_vector();
      BOOST_FOREACH(real value, result_real_vector) ret.append(value);
      return ret;
    }

    python::list Observable::compute_int_vector_python() {
      python::list ret;
      compute_int_vector();
      BOOST_FOREACH(int value, result_int_vector) ret.append(value);
      return ret;
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Observable::registerPython() {
      using namespace espressopp::python;
    
      class_< Observable, boost::noncopyable >("analysis_Observable", no_init)
		.def("compute", &Observable::compute)
		.def("compute_real", &Observable::compute_real)
		.def("compute_int", &Observable::compute_int)
		.def("compute_real_vector_python", &Observable::compute_real_vector_python)
		.def("compute_int_vector_python",  &Observable::compute_int_vector_python)
		.def("getResultType", &Observable::getResultType)
      ;
    }
  }
}
