#include "Timer.hpp"

#ifdef HAVE_SYS_RESOURCE_H

#include <sys/resource.h>

using namespace espresso::esutil;

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
