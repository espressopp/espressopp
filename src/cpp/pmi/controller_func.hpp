#ifndef _PMI_CONTROLLER_FUNC_HPP
#define _PMI_CONTROLLER_FUNC_HPP

#include "pmi/types.hpp"

#ifdef CONTROLLER
namespace pmi {
  bool isWorkersActive();
  void endWorkers();
}

#endif /* CONTROLLER */
#endif /* _PMI_CONTROLLER_HPP */
