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
#include "AnalysisBase.hpp"

//#include "ParticleAccess.hpp"

namespace espressopp {
  namespace analysis {

    LOG4ESPP_LOGGER(AnalysisBase::logger, "AnalysisBase");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    AnalysisBase::registerPython() {
      using namespace espressopp::python;
    
      //class_< AnalysisBase, boost::noncopyable >("analysis_AnalysisBase", no_init)
      class_< AnalysisBase, bases<ParticleAccess>, boost::noncopyable >("analysis_AnalysisBase", no_init)
          .def("performMeasurement", pure_virtual(&AnalysisBase::performMeasurement))
	  .def("reset", pure_virtual(&AnalysisBase::reset))
	  .def("compute", pure_virtual(&AnalysisBase::compute))
	  .def("getAverageValue", pure_virtual(&AnalysisBase::getAverageValue))
	  .def("getNumberOfMeasurements", pure_virtual(&AnalysisBase::getNumberOfMeasurements))
      ;
    }
  }
}
