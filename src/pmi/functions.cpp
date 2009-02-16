#include "pmi/functions.hpp"
#include "pmi/worker_internal.hpp"

#include <sstream>

using namespace std;

namespace pmi {
  unsigned int getControllerMPIRank() {
    return CONTROLLER_ID;
  }

  WorkerIdType getWorkerId() {
    return pmi::transmit::getWorkerId();
  }

  bool 
  isController() {
    return getWorkerId() == CONTROLLER_ID;
  }
  
  bool 
  isWorker() {
    return !isController();
  }
  
  // pretty print the worker id
  string
  printWorkerId() {
    if (isController())
      return "Controller ";
    else {
      ostringstream ost;
      ost << "Worker " << getWorkerId() << " ";
      return ost.str();
    }
  }
  
  bool 
  mainLoop() {
    if (isController()) {
      LOG4ESPP_INFO(logger, "Controller starts and exits main loop.");
      return false;
    }

    LOG4ESPP_INFO(logger, "Worker " << getWorkerId() << " starts main loop.");
    
    // The main loop
    while (transmit::handleNext()) {};

    LOG4ESPP_INFO(logger, "Worker " << getWorkerId() << " finished.");
    return true;
  }

  bool workersActive = true;

  bool 
  isWorkersActive() { return workersActive; }
  
  void 
  endWorkers() {
    if (isWorker())
      PMI_THROW_USER_ERROR(printWorkerId() << "tries to end workers.");
    LOG4ESPP_INFO(logger, "Controller ends all workers.");
    transmit::endWorkers();
    workersActive = false;
#ifndef PMI_OPTIMIZE
    transmit::gatherStatus();
#endif
  }
}



