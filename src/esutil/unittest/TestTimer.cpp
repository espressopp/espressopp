#include <acconfig.hpp>
#define BOOST_TEST_MODULE timing
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../Timer.hpp"

using namespace esutil;

int doSomething()
{
    const int size = 1000;
    int cnt = 0;
    for (int x = 1; x < size; ++x) {
        for (int y = 1; y < size; ++y) {
            for (int z = 1; z < size; ++z) {
                if (x*x + y*y == z*z) {
                    cnt++;
                }
            }
        }
    }
    return cnt;
}

BOOST_AUTO_TEST_CASE(timing_test) {
    UserTimer utimer;
    WallTimer wtimer;

    doSomething();
    
    BOOST_TEST_MESSAGE("user time " << utimer);
    BOOST_TEST_MESSAGE("wall time " << wtimer);

    // basically, everything between 0 and 20 secs is reasonable,
    // just not zero, if we have reasonable precision
    BOOST_WARN_MESSAGE(utimer.getElapsedTime() > 0, "user time is zero, probably not implemented for your platform. Timing will not work.");
    BOOST_WARN_MESSAGE(utimer.getElapsedTime() < 20, "user time is exceedingly large, either your platform is slow or timing broken.");
    BOOST_WARN_MESSAGE(wtimer.getElapsedTime() > 0, "wall time is zero, probably not implemented for your platform. Timing will not work.");
    BOOST_WARN_MESSAGE(wtimer.getElapsedTime() < 20, "wall time is exceedingly large, either your platform is slow or timing broken.");
}

/*
  Local Variables:
  compile-command: "mpic++ -Wall -g -I. -I../.. \
  -I/home/axel/software/include/boost-1_36 \
  -L/home/axel/software/lib Timer.cpp \
  ../Timer.cpp -o timing \
  -lboost_unit_test_framework-gcc42-mt-1_36 -lboost_mpi-gcc42-mt-1_36 -lboost_serialization-gcc42-mt-1_36 && ./timing --log_level=message"
  End:
*/
