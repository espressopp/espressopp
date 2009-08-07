#ifndef _THERMOSTAT_LANGEVIN_HPP
#define _THERMOSTAT_LANGEVIN_HPP

#include <boost/signals2.hpp>
#include <boost/random.hpp>

#include "logging.hpp"
#include "types.hpp"
#include "thermostat/Thermostat.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "Real3D.hpp"

namespace espresso {
  namespace thermostat {

    /** A Langevin thermostat models the embedding of particles into a
       viscouls liquid of constant temperature using kicks and drag.

       The Langevin thermostat is connected to a Velocity-Verlet 
       integrator (\sa espresso::integrator::VelocityVerlet).

    */

    class Langevin: 
      public Thermostat, 
      public enable_shared_from_this< Langevin >
    {

    public:
      typedef shared_ptr< Langevin > SelfPtr;

    private:
  
      static LOG4ESPP_DECL_LOGGER(theLogger);

      real gamma;

      // random number generator for normal distribution

      boost::minstd_rand0 linearCongruential;

      boost::normal_distribution<double> normalDist;

      boost::variate_generator<boost::minstd_rand0&, 
			       boost::normal_distribution<double> > gauss;

      /** Method of the thermostat to modify r and v before force computation. */
      virtual void thermalizeA(const integrator::VelocityVerlet&);

      /** Method of the thermostat to modify r and v after force computation. */
      virtual void thermalizeB(const integrator::VelocityVerlet&);

    public:

      /** Method to register the Python bindings. */
      static void registerPython();

      /** Constructor with the particle set, temperature, and friction coefficient. */
      Langevin(real temperature, real gamma);

      /** The thermostat will be connected to a VelocityVerlet integrator so that it
          can be called by this integrator. 

          \param thisSharedPtr must be the shared pointer to this object
          \param integrator is the Velocity-Verlet integrator to which this thermostat will be connected.

          This thermostat will not keep a reference to the integrator but to its connections.

      */

      void connect(integrator::VelocityVerlet::SelfPtr integrator);

      /** Disconnect this thermostat object from its integrator. */
      void disconnect(integrator::VelocityVerlet::SelfPtr integrator);

      /** Getter for gamma. */
      real getGamma() const;

      /** Setter for gamma. */
      void setGamma(real gamma);

      virtual ~Langevin();
    };
  }
}

#endif
