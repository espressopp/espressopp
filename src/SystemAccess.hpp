#ifndef _SYSTEM_ACCESS_HPP
#define _SYSTEM_ACCESS_HPP

#include "types.hpp"
#include "System.hpp"

namespace espresso {

  /** Common base class for all classes that need access to the system.

      This class encapsulates the access to a system via weak pointers and
      takes care of error handling (NULL system, system expired)

  */

  class SystemAccess {

   public:

    /** Constructor */

    SystemAccess(shared_ptr<System> system);

    /** Get shared pointer to the system */

    shared_ptr<System> getSystem() const;

  private:

    weak_ptr< System > mySystem;

  };

  /**************************************************************************
  *  Implementation                                                         *
  **************************************************************************/

  inline SystemAccess::SystemAccess(shared_ptr<System> system) {
    if (!system) {
       throw std::runtime_error("NULL system");
    } 

    // This is a workaround currently needed for Boost-Python
    // This does not work: mySystem = system

    if (!system->getShared()) {
       throw std::runtime_error("INTERNAL error: no shared pointer for system");
    }

    mySystem = system->getShared();
  }

  inline shared_ptr<System> SystemAccess::getSystem() const {
 
    if (mySystem.expired()) {
       throw std::runtime_error("expired system");
    }

    return mySystem.lock();
  }
}

#endif
