/*
 Copyright (C) 2012-2016 Max Planck Institute for Polymer Research
 Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

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

#include "LBOutput.hpp"

namespace espressopp {
   namespace analysis {
      //    LOG4ESPP_LOGGER(LBOutput::theLogger, "LBOutput");
      /* Abstract class. No need to define anything here, except of python registration */
      void LBOutput::registerPython() {
         using namespace espressopp::python;

         class_<LBOutput, bases< AnalysisBase >, boost::noncopyable >
         ("analysis_LBOutput", no_init)

         .def("writeOutput", pure_virtual(&LBOutput::writeOutput))
         ;
      }
   }
}
