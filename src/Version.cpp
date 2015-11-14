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
#include "Version.hpp"
#include <sstream>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
# define VT_ON()
# define VT_OFF()
#endif

namespace espressopp {

  Version::Version() {
	  name          = "ESPResSo++";
      major         = MAJORVERSION;
	  minor         = MINORVERSION;
	  patchlevel    = PATCHLEVEL;
	  gitrevision    = gitversion;
      boostversion  = BOOST_LIB_VERSION;
	  date          = __DATE__;
	  time          = __TIME__;
  }

  std::string Version::info() {
	  std::stringstream ss;
	  ss << name << " v"  << major << "." << minor;
	  ss << " patchlevel " << patchlevel;
	  if (gitrevision != "")
	    ss << ", Git revision: " << gitrevision;
      ss << ", Boost Version: " << boostversion;
	  ss << ", compiled on " << date << ", " << time;
#ifdef VTRACE
	  ss << ", VampirTrace mode";
#endif
	  return ss.str();
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  Version::registerPython() {
    using namespace espressopp::python;

    class_< Version >
      ("Version", init<>())
      .def_readonly("major", &Version::major)
      .def_readonly("minor", &Version::minor)
      .def_readonly("gitrevision", &Version::gitrevision)
      .def_readonly("boostversion", &Version::boostversion)
      .def_readonly("patchlevel", &Version::patchlevel)
      .def_readonly("date", &Version::date)
      .def_readonly("time", &Version::time)
      .def_readonly("name", &Version::name)
      .def("info", &Version::info);
  }
}
