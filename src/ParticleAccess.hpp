#ifndef _PARTICLEACCESS_HPP
#define	_PARTICLEACCESS_HPP
#include "types.hpp"
#include "python.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  class ParticleAccess : public SystemAccess {
  public:
    ParticleAccess(shared_ptr< System > system) : SystemAccess(system) {}
    virtual ~ParticleAccess() {}

    virtual void perform_action() = 0;
    
    static void registerPython();
  };

}
#endif	/* PARTICLEACCESS_HPP */