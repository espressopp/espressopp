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

// ESPP_CLASS
#ifndef _VERSION_HPP
#define _VERSION_HPP
#include <string>
#include "boost/version.hpp"

#define MAJORVERSION 2
#define MINORVERSION 0
#define PATCHLEVEL   0
#include "gitversion.hpp"

namespace espressopp {

  class Version {
  public:
    Version();
    std::string info();
    static void registerPython();

  private:
    int major;
    int minor;
    int patchlevel;
    std::string name;
    std::string gitrevision;
    std::string boostversion;
    std::string date;
    std::string time;
  };
}
#endif
