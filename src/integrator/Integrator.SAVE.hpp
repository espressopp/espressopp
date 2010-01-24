#ifndef _INTEGRATOR_HPP
#define _INTEGRATOR_HPP

//////////////////////////////////////////////////
// Mock classes, expected interfaces
//////////////////////////////////////////////////
#include "types.hpp"

class Error {};

class ErrorHandler {
  void reportError(Error error) {
    // collect errors
  }
  bool checkErrors() {
    // gather_all error status from nodes
    // return whether an error was found on any node
  }
  list< Error > getErrors() {
    // get the collected errors (only on the master node)
  }
};

class System {
  ErrorHandler::SelfPtr errorHandler;

public:
  typedef boost::shared_ptr< System > SelfPtr;

  ErrorHandler::SelfPtr getErrorHandler() {
    return errorHandler;
  }
};

// Signals:
// * startIntegration
// Dynamic registry:
// * 

// * Thermostat has to connect to signal startIntegration
//   * uses flags: reinit_thermo and recalc_forces

//////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////
#include <boost/signals2.hpp>

#include "logging.hpp"

class VelocityVerlet {
  System::SelfPtr system;

  static LOG4ESPP_DECL_LOGGER(logger);

public:
  boost::signals2::signal1< void > startIntegration;
  
  /**
     Check whether the integrator has valid values for its parameters.
  */
  bool checkSanity() {
    
  }

  void integrate(int n_steps) {
    checkSanity();
    system->getErrorHandler->checkErrors();

    // call signal
    startIntegration();

    // TODO: Should go into verlet list implementation
    // /* Verlet list criterion */
    // skin2 = SQR(0.5 * skin);

    LOG4ESPP_INFO(logger, "Integrating " << n_steps << " steps.");

    
  }
};

#endif
