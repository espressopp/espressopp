#ifndef _PMI_TYPES_HPP
#define _PMI_TYPES_HPP

#include "acconfig.hpp"
#include "pmi/log4espp.hpp"
#include <limits>
#include <string>
#include <exception>

#if !defined(CONTROLLER) && !defined(WORKER)
#define CONTROLLER 1
#define WORKER 1
#endif

extern LOG4ESPP_DECL_LOGGER(logger);

namespace pmi {
  // class required to call class initializers
  class Nothing {};

  typedef unsigned int IdType;
  const static IdType NOT_ASSOCIATED = 
    std::numeric_limits<IdType>::max();
  const static IdType NOT_DEFINED = 
    std::numeric_limits<IdType>::max();

  const static std::string NOT_REGISTERED =
    "NOT_REGISTERED";

  typedef unsigned int WorkerIdType;
  const static int CONTROLLER_ID = 0;
}

#endif
