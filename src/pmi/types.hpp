#ifndef _PMI_TYPES_HPP
#define _PMI_TYPES_HPP

#include "acconfig.hpp"
#include <log4espp.hpp>
#include <limits>
#include <string>
#include <exception>

namespace pmi {
  extern LOG4ESPP_DECL_LOGGER(logger);

  typedef unsigned int IdType;
  const static IdType NOT_ASSOCIATED = 
    std::numeric_limits<IdType>::max();
  const static IdType NOT_DEFINED = 
    std::numeric_limits<IdType>::max();

  const static std::string NOT_REGISTERED =
    "NOT_REGISTERED";

  typedef unsigned int WorkerIdType;
  const static WorkerIdType CONTROLLER_ID = 0;
}

#endif
