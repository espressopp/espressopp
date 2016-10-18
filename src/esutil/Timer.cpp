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

#include "Timer.hpp"

namespace espressopp {
namespace esutil {

#ifdef CMAKE_HEADERS

#include <sys/resource.h>

float UserTimer::getCurrentTime() const {
  struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec + 1.e-6*rus.ru_utime.tv_usec;
}

#else

// we do not have getrusage
float UserTimer::getCurrentTime() const { return 0; }

#endif

float WallTimer::getCurrentTime() const { return timer.elapsed(); }

}  // end namespace esutil
}  // end namespace espressopp
