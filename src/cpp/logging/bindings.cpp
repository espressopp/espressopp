#include "logging.hpp"

// TODO: fix log4espp to make it possible to move
// this to the logging namespace
LOG4ESPP_DEFINITION();

namespace logging {

  void initLogging()
  {
    LOG4ESPP_CONFIGURE();
  }
}
