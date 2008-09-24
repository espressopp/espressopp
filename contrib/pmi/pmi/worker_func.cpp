#include "pmi/worker_func.hpp"
#include "pmi/worker_internal.hpp"
#include "pmi/transmit.hpp"

#ifdef WORKER

using namespace std;

namespace pmi {
  // These store the callers between registration and association
  // construct-on-first-use idiom
  map<string, ConstructorCallerType> &constructorCallersByName() {
    static map<string, ConstructorCallerType> callers;
    return callers;
  }

  // construct-on-first-use idiom
  map<string, MethodCallerType> &methodCallersByName() {
    static map<string, MethodCallerType> callers;
    return callers;
  }

  // construct-on-first-use idiom
  map<string, DestructorCallerType> &destructorCallersByName() {
    static map<string, DestructorCallerType> callers;
    return callers;
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
}
#endif

