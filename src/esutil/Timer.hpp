#ifndef _ESUTIL_TIMER_HPP
#define _ESUTIL_TIMER_HPP

#include "acconfig.hpp"

#include <ostream>
#include <boost/mpi/timer.hpp>

namespace espresso {
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
