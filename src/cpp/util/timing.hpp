#ifndef _UTIL_TIMING_HPP
#define _UTIL_TIMING_HPP
#include "acconfig.hpp"
#include "types.hpp"

namespace util {
    /** get the current user time in secs. The starting time is
	arbitrary, therefore use this function only to calculate time
	spans (differences). */
    real userSecs();
    
    /** get the current wall time in secs. The starting date is
	arbitrary, therefore use this function only to calculate time
	spans (differences). */
    real wallSecs();
}
#endif
