// ESPP_CLASS
#ifndef _ADRESS_HPP
#define _ADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "integrator/VelocityVerlet.hpp"

#include "boost/signals2.hpp"


namespace espresso {

  class Adress : public SystemAccess {

  public:

    Adress(shared_ptr<System> system,
            shared_ptr<integrator::VelocityVerlet> _integrator);

    ~Adress();


    /** Register this class so it can be used from Python. */
    static void registerPython();

  private:


    boost::signals2::connection _initForces, _integrate1, _integrate2;

    shared_ptr<integrator::VelocityVerlet> integrator; // this is needed for signal connection


    void integrate1(real&);
    void initForces();
    void integrate2();


    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
