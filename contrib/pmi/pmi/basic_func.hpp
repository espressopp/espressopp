#ifndef _PMI_BASIC_FUNC_HPP
#define _PMI_BASIC_FUNC_HPP

#include "pmi/types.hpp"
#include "pmi/transmit.hpp"

namespace pmi {
  // test whether the executing process is the controller
  bool isController();
  // test whether the executing process is the worker
  bool isWorker();
  // pretty print the worker id
  std::string printWorkerId();
}
#endif
