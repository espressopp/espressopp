#include "python.hpp"
#include "Version.hpp"
#include <sstream>

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
# define VT_ON()
# define VT_OFF()
#endif

namespace espresso {

  Version::Version() {
	  name          = "ESPResSo++";
      major         = MAJORVERSION;
	  minor         = MINORVERSION;
	  patchlevel    = PATCHLEVEL;
	  hgrevision    = hgversion;
      boostversion  = BOOST_LIB_VERSION;
	  date          = __DATE__;
	  time          = __TIME__;
  }

  std::string Version::info() {
	  std::stringstream ss;
	  ss << name << " v"  << major << "." << minor;
	  ss << " patchlevel " << patchlevel;
	  if (hgrevision != "")
	    ss << ", Mercurial(hg) revision: " << hgrevision;
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
    using namespace espresso::python;

    class_< Version >
      ("Version", init<>())
      .def_readonly("major", &Version::major)
      .def_readonly("minor", &Version::minor)
      .def_readonly("hgrevision", &Version::hgrevision)
      .def_readonly("boostversion", &Version::boostversion)
      .def_readonly("patchlevel", &Version::patchlevel)
      .def_readonly("date", &Version::date)
      .def_readonly("time", &Version::time)
      .def_readonly("name", &Version::name)
      .def("info", &Version::info);
  }
}
