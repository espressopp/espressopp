#include "timing.hpp"

#ifdef HAVE_SYS_RESOURCE_H

#include <sys/resource.h>

namespace util {
    real userSecs() {
        struct rusage rus;
        getrusage(RUSAGE_SELF, &rus);
        return rus.ru_utime.tv_sec + 1.e-6*rus.ru_utime.tv_usec;
    }
}

#else

namespace util {
    // we do not have getrusage
    real userSecs() { return 0; }
}
#endif

#ifdef HAVE_SYS_TIME_H

#include <sys/time.h>

namespace util {
    real wallSecs() {
        struct timeval tp;
        struct timezone tzp;
        gettimeofday (&tp, &tzp);
        return tp.tv_sec + 1.e-6*tp.tv_usec;
    }
}

#elif defined(HAVE_TIME_H)

#include <time.h>

namespace util {
    real wallSecs() {
        return time(0);
    }
}

#else

namespace util {
    // we do not have gettimeofday
    real wallSecs() { return 0; }
}

#endif
