
/** ARGERROR specifies illegal argument for any routine, occurs on all processors at same time

    An exception is thrown that will be caught at Python script level.
*/

#include <stdexcept>
#include "logging.hpp"

#define ARGERROR(logger,msg) { std::ostringstream omsg; omsg << msg; \
                     LOG4ESPP_ERROR(logger, msg); throw std::invalid_argument(omsg.str()); }

