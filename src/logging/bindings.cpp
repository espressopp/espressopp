#include "logging.hpp"

namespace logging {

  void initLogging()
  {
    LOG4ESPP_CONFIGURE();
  }

  void finalizeLogging()
  {
  }
}
