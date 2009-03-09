#ifndef _ESUTIL_TIMER_HPP
#define _ESUTIL_TIMER_HPP

#include "acconfig.hpp"

#include <ostream>
#ifdef HAVE_BOOST_MPI
#include <boost/mpi/timer.hpp>
#endif

#include "estypes.hpp"

namespace espresso {
  namespace esutil {
    /** simple timer for obtaining typically microsecond precision
        timings.  The time is measured in seconds from the point of
        construction or whenever reset() is called. The precision
        depends on the available time measuring functions.
    */
    class Timer {
    protected:
      real currentTime;
      virtual real getCurrentTime() const = 0;
    public:
      virtual ~Timer() {}
      /// reset the starting time
      void reset() { currentTime = getCurrentTime(); }
      /// get the time that elapsed since the last reset
      real getElapsedTime() const { return getCurrentTime() - currentTime; }
    };

    /// when printing give the current elapsed time
    inline std::ostream &operator<<(std::ostream &os, const Timer &timer) {
      os << timer.getElapsedTime() << "s";
      return os;
    }

    /// timer measuring the user time.
    class UserTimer: public Timer {
      virtual real getCurrentTime() const;
    public:
      UserTimer() { reset(); }
    };

#ifdef HAVE_BOOST_MPI
 
    /// timer measuring the wall time.
    class WallTimer: public Timer {
      boost::mpi::timer timer;

      virtual real getCurrentTime() const;
    };

#else

    /// timer measuring the wall time.
    class WallTimer: public Timer {
      virtual real getCurrentTime() const;
    public:
      WallTimer() { reset(); }
    };

#endif
  }
}
#endif
