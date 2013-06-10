#ifndef _SYSTEM_ACCESS_HPP
#define _SYSTEM_ACCESS_HPP

#include "System.hpp"
#include "types.hpp"

namespace espresso {

  /** Common base class for all classes that need access to the system.

      This class encapsulates the access to a system via weak pointers and
      takes care of error handling (NULL system, system expired)
  */

  class SystemAccess {
  public:
    /** Constructor */
    SystemAccess(shared_ptr<System> system);
    SystemAccess(shared_ptr<System> system1, shared_ptr<System> system2);

    /** Get shared pointer to the system */
    shared_ptr<System> getSystem() const;
    shared_ptr<System> getSystem1() const;
    shared_ptr<System> getSystem2() const;

    /** Get a system reference */
    System& getSystemRef() const { return *(getSystem().get()); }
    System& getSystemRef1() const { return *(getSystem1().get()); }
    System& getSystemRef2() const { return *(getSystem2().get()); }
  private:
    /** Weak pointer to the system guarantees that system will be freed */
    weak_ptr< System > mySystem;
    weak_ptr< System > mySystem1;
    weak_ptr< System > mySystem2;
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

  inline SystemAccess::SystemAccess(shared_ptr<System> system1, shared_ptr<System> system2) {
	if (!system1) {
	  throw std::runtime_error("NULL system1");
	}
	if (!system2) {
	  throw std::runtime_error("NULL system2");
    }
    // This is a workaround currently needed for Boost-Python
    // This does not work: mySystem = system
    if (!system1->getShared()) {
       throw std::runtime_error("INTERNAL error: no shared pointer for system1");
    }
    if (!system2->getShared()) {
       throw std::runtime_error("INTERNAL error: no shared pointer for system2");
    }

    mySystem1 = system1->getShared();
    mySystem2 = system2->getShared();
  }


  inline shared_ptr<System> SystemAccess::getSystem() const {
 
    if (mySystem.expired()) {
       throw std::runtime_error("expired system");
    }

    return mySystem.lock();
  }

  inline shared_ptr<System> SystemAccess::getSystem1() const {

    if (mySystem1.expired()) {
       throw std::runtime_error("expired system1");
    }

    return mySystem1.lock();
  }

  inline shared_ptr<System> SystemAccess::getSystem2() const {

    if (mySystem2.expired()) {
       throw std::runtime_error("expired system2");
    }

    return mySystem2.lock();
  }

}

#endif
