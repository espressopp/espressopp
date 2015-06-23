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

#include <string>
#include <sstream>
#include "python.hpp"
#include "Version.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
# define VT_ON()
# define VT_OFF()
#endif

namespace espressopp {

  Version::Version() {
    name_          = "ESPResSo++";
    major_         = MAJORVERSION;
    minor_         = MINORVERSION;
    patchlevel_    = PATCHLEVEL;
    gitrevision_    = gitversion;
    boostversion_  = BOOST_LIB_VERSION;
    date_          = __DATE__;
    time_          = __TIME__;
  }

  std::string Version::info() {
    std::stringstream ss;
    ss << name_;
    ss << version();
    ss << " Boost Version: " << boostversion_;
    ss << ", compiled on " << date_ << ", " << time_;
#ifdef VTRACE
    ss << ", VampirTrace mode";
#endif
    return ss.str();
  }

  std::string Version::version() {
    std::stringstream ss;
    ss << major_ << "." << minor_;
    ss << " (patchlevel " << patchlevel_;
    if (gitrevision_ != "")
      ss << ", revision: " << gitrevision_;
    ss << ")";
    return ss.str();
  }

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  Version::registerPython() {
    using namespace espressopp::python;  //NOLINT

    class_< Version >
      ("Version", init<>())
      .def_readonly("major", &Version::major_)
      .def_readonly("minor", &Version::minor_)
      .def_readonly("gitrevision", &Version::gitrevision_)
      .def_readonly("boostversion", &Version::boostversion_)
      .def_readonly("patchlevel", &Version::patchlevel_)
      .def_readonly("date", &Version::date_)
      .def_readonly("time", &Version::time_)
      .def_readonly("name", &Version::name_)
      .def("info", &Version::info);
  }
}  // end namespace espressopp
