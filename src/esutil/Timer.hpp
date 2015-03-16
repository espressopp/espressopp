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
#ifndef _ESUTIL_TIMER_HPP
#define _ESUTIL_TIMER_HPP

#include "acconfig.hpp"

#include <ostream>
#include <boost/mpi/timer.hpp>

namespace espressopp {
  namespace esutil {
    /** simple timer for obtaining typically microsecond precision
        timings.  The time is measured in seconds from the point of
        construction or whenever reset() is called. The precision
        depends on the available time measuring functions.
    */
    class Timer {
    protected:
      float currentTime;
      virtual float getCurrentTime() const = 0;
    public:
      virtual ~Timer() {}
      /// reset the starting time
      void reset() { currentTime = getCurrentTime(); }
      /// get the time that elapsed since the last reset
      float getElapsedTime() const { return getCurrentTime() - currentTime; }
    };

    /// when printing give the current elapsed time
    inline std::ostream &operator<<(std::ostream &os, const Timer &timer) {
      os << timer.getElapsedTime() << "s";
      return os;
    }

    /// timer measuring the user time.
    class UserTimer: public Timer {
      virtual float getCurrentTime() const;
    public:
      UserTimer() { reset(); }
    };

    /// timer measuring the wall time.
    class WallTimer: public Timer {
      boost::mpi::timer timer;

      virtual float getCurrentTime() const;
    };
  }
}
#endif
